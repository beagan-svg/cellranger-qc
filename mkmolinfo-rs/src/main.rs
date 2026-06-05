//! mkmolinfo — Rust implementation of the `mkmolinfo` utility.
//!
//! The CLI, input/output filenames, HDF5 dataset names, 2-bit barcode/UMI
//! encoder, per-molecule metric layout, and per-read classification branches are
//! defined for this utility. Counter datasets are u32 and barcode
//! encodings are u64.
//!
//! Purpose:
//!   Takes a path to a Cell Ranger output directory and creates
//!   `molecule_info_new.h5`.
//!
//! High-level flow:
//!   1. Locate `possorted_genome_bam.bam` and `molecule_info.h5` inside the
//!      given Cell Ranger output directory.
//!   2. Copy the input h5 to the output path (default `molecule_info_new.h5`)
//!      and open the copy read-write. Genome and metadata groups carry over.
//!   3. Read `/barcode_idx`, `/umi` and `/barcodes` from the copied file
//!      (three `Container::read_1d` calls). These define the molecule list.
//!   4. Verify every barcode is from a single GEM group (suffix "-1"); the command
//!      exits otherwise ("Barcode is from multiple GEM groups...").
//!   5. Stream the BAM. Per record the tags are read in this exact order:
//!      CB, UB, CR, UR, TX. A read is joined to a molecule by
//!      (encode(CB), encode(UB)); only molecules present in the file are kept.
//!   6. Write exactly SEVEN datasets back into the copied file (dataset names):
//!      reads, barcode_corrected_reads, conf_mapped, nonconf_mapped_reads,
//!      umi_corrected_reads, unmapped_reads  (six u32 arrays via `write_array`),
//!      and `barcode` (one u64 array, written directly).
//!      NOTE: this command does not write gene/genome/gem_group/umi — those already
//!      live in the copied file.
//!
//! Per-read classification:
//!   entry.reads                   += 1   // [rbx+0x00] every joined read
//!   if CR present && CR != CB_seq:
//!       entry.barcode_corrected   += 1   // [rbx+0x04]
//!   if UR present && encode(UR) != encode(UB):
//!       entry.umi_corrected       += 1   // [rbx+0x10]
//!   if tid == -1:
//!       entry.unmapped            += 1   // [rbx+0x14]
//!   else if mapq != 255:
//!       entry.nonconf_mapped      += 1   // [rbx+0x0c]
//!   else if TX tag present:              // tid != -1 && mapq == 255
//!       entry.conf_mapped         += 1   // [rbx+0x08]
//!   // a uniquely-mapped read (mapq 255) lacking a TX tag counts in `reads` only.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use clap::Parser;
use ndarray::Array1;
use rust_htslib::bam::{self, Read};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Default output filename.
const DEFAULT_OUTPUT: &str = "molecule_info_new.h5";
/// Input BAM filename searched for in the output dir.
const BAM_FILENAME: &str = "possorted_genome_bam.bam";
/// Input molecule info to copy from (searched for in the dir).
const MOLECULE_INFO_FILENAME: &str = "molecule_info.h5";

// BAM aux tags read per record.
const TAG_CB: &[u8] = b"CB"; // corrected cell barcode (carries the "-1" GEM suffix)
const TAG_UB: &[u8] = b"UB"; // corrected UMI
const TAG_CR: &[u8] = b"CR"; // raw (uncorrected) cell barcode
const TAG_UR: &[u8] = b"UR"; // raw (uncorrected) UMI
const TAG_TX: &[u8] = b"TX"; // transcriptome-compatible alignment (gates conf_mapped)

/// In 10x BAMs a uniquely / confidently mapped read carries MAPQ 255.
const CONF_MAPPED_MAPQ: u8 = 255;

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

#[derive(Parser, Debug)]
#[command(
    name = "mkmolinfo",
    version = "1.0",
    about = "Takes a path to a Cell Ranger output directory and creates molecule_info_new.h5."
)]
struct Cli {
    /// Specify the output directory produced by Cell Ranger
    #[arg(value_name = "OUT_DIRECTORY")]
    out_directory: PathBuf,

    /// Specify an output file name (default molecule_info_new.h5)
    #[arg(short = 'o', long = "output", default_value = DEFAULT_OUTPUT)]
    output: PathBuf,

    /// Print progress as program parses through BAM file.
    #[arg(short = 'v')]
    verbose: bool,

    /// Diagnostic: run the join (steps 1-5), report how many
    /// molecules from the file were found in the BAM, then exit without writing.
    #[arg(long = "self-check")]
    self_check: bool,
}

// ---------------------------------------------------------------------------
// Barcode/UMI 2-bit encoder.
// The asm subtracts 'A' (0x41) and uses a jump table, shifting left 2 bits per
// base: rax = (rax << 2) | code, with A=0,C=1,G=2,T=3 (MSB-first). The same
// routine encodes both cell barcodes and UMIs (called twice per record on
// CB and UB). A "-1" GEM-group suffix on CB is dropped before encoding.
// ---------------------------------------------------------------------------

fn barcode_str_to_u64(seq: &str) -> Result<u64> {
    let seq = seq.split('-').next().unwrap_or(seq);
    let mut acc: u64 = 0;
    for b in seq.bytes() {
        let code: u64 = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            other => bail!(
                "Could not read in barcode strings: invalid base {:?}",
                other as char
            ),
        };
        acc = (acc << 2) | code;
    }
    Ok(acc)
}

/// The barcode sequence portion (everything before the "-<gem>" suffix).
fn barcode_seq(barcode: &str) -> &str {
    barcode.split('-').next().unwrap_or(barcode)
}

/// GEM-group suffix of a barcode string ("ACGT...-1" -> 1).
fn gem_group_of(barcode: &str) -> Option<u16> {
    barcode
        .rsplit('-')
        .next()
        .and_then(|s| s.parse::<u16>().ok())
}

// ---------------------------------------------------------------------------
// Per-molecule accumulator — 6 × u32, matching the metric layout.
// Field order here keeps the field order explicit for clarity.
// ---------------------------------------------------------------------------

#[derive(Default, Clone)]
struct MoleculeMetrics {
    reads: u32,                   // [rbx+0x00]
    barcode_corrected_reads: u32, // [rbx+0x04]
    conf_mapped: u32,             // [rbx+0x08]
    nonconf_mapped_reads: u32,    // [rbx+0x0c]
    umi_corrected_reads: u32,     // [rbx+0x10]
    unmapped_reads: u32,          // [rbx+0x14]
}

/// Molecule identity used to join BAM reads to molecule_info rows:
/// (encoded cell barcode, encoded UMI).
type MolKey = (u64, u64);

// ---------------------------------------------------------------------------
// HDF5 helpers
// ---------------------------------------------------------------------------

/// Generic dataset writer — implementation of `mkmolinfo::write_array`
/// (`DatasetBuilder::create` followed by `Container::write`). This helper has a
/// single monomorphization (u32), used for all six count datasets.
fn write_array<T>(file: &hdf5::File, name: &str, data: &Array1<T>) -> Result<()>
where
    T: hdf5::H5Type,
{
    let ds = file
        .new_dataset::<T>()
        .shape([data.len()])
        .create(name)
        .with_context(|| format!("failed to create dataset {name}"))?;
    ds.write(data)
        .with_context(|| format!("failed to write dataset {name}"))?;
    Ok(())
}

/// Read a 1-D integer dataset as u64 (handles the widths molecule_info uses).
fn read_u64_dataset(file: &hdf5::File, name: &str, missing_msg: &str) -> Result<Vec<u64>> {
    let ds = file.dataset(name).map_err(|_| anyhow!("{missing_msg}"))?;
    if let Ok(v) = ds.read_1d::<u64>() {
        return Ok(v.to_vec());
    }
    if let Ok(v) = ds.read_1d::<u32>() {
        return Ok(v.into_iter().map(|x| x as u64).collect());
    }
    if let Ok(v) = ds.read_1d::<i64>() {
        return Ok(v.into_iter().map(|x| x as u64).collect());
    }
    let v = ds
        .read_1d::<i32>()
        .with_context(|| format!("Could not read in {name}"))?;
    Ok(v.into_iter().map(|x| x as u64).collect())
}

/// Read the `/barcodes` string table.
fn read_barcode_strings(file: &hdf5::File) -> Result<Vec<String>> {
    let ds = file
        .dataset("/barcodes")
        .map_err(|_| anyhow!("H5 file did not contain barcodes"))?;
    if let Ok(v) = ds.read_1d::<hdf5::types::FixedAscii<32>>() {
        return Ok(v.iter().map(|s| s.as_str().to_string()).collect());
    }
    let v = ds
        .read_1d::<hdf5::types::VarLenAscii>()
        .context("Could not read in barcode strings")?;
    Ok(v.iter().map(|s| s.as_str().to_string()).collect())
}

// ---------------------------------------------------------------------------
// BAM streaming + classification
// ---------------------------------------------------------------------------

/// Pull an aux tag as a string, if present and string-typed.
fn aux_string(rec: &bam::Record, tag: &[u8]) -> Option<String> {
    match rec.aux(tag) {
        Ok(bam::record::Aux::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

/// True if the record carries the given aux tag at all (any type).
fn has_aux(rec: &bam::Record, tag: &[u8]) -> bool {
    rec.aux(tag).is_ok()
}

/// Parse the BAM, accumulating per-molecule metrics keyed by (barcode, umi).
/// `valid_keys` restricts work to molecules present in molecule_info.h5.
fn parse_bam(
    bam_path: &Path,
    valid_keys: &HashSet<MolKey>,
    verbose: bool,
) -> Result<HashMap<MolKey, MoleculeMetrics>> {
    let mut reader = bam::Reader::from_path(bam_path)
        .with_context(|| format!("Could not open BAM {}", bam_path.display()))?;

    let mut metrics: HashMap<MolKey, MoleculeMetrics> = HashMap::new();
    let mut record = bam::Record::new();
    let mut n: u64 = 0;

    while let Some(res) = reader.read(&mut record) {
        res.context("error reading BAM record")?;
        n += 1;
        if verbose && n % 5_000_000 == 0 {
            eprintln!("Parsed {n} BAM records...");
        }

        // Tags are read in a fixed order: CB, UB, (CR, UR, TX).
        let cb = match aux_string(&record, TAG_CB) {
            Some(s) => s,
            None => continue, // no corrected barcode -> skip record
        };
        // The command checks that the CB's last byte is '1' (GEM group 1).
        // Any other GEM group is a fatal error — it aborts the
        // whole run rather than skipping the read.
        if cb.as_bytes().last() != Some(&b'1') {
            bail!(
                "Error: Barcode is from multiple GEM groups, only one is supported. \
                 Please email the author to enable this. BC: {cb}"
            );
        }
        let ub = match aux_string(&record, TAG_UB) {
            Some(s) => s,
            None => continue,
        };

        let bc = match barcode_str_to_u64(&cb) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let umi = match barcode_str_to_u64(&ub) {
            Ok(v) => v,
            Err(_) => continue,
        };

        let key: MolKey = (bc, umi);
        if !valid_keys.contains(&key) {
            continue;
        }

        let m = metrics.entry(key).or_default();
        m.reads += 1;

        // Barcode correction: raw CR vs the CB *sequence* (suffix stripped).
        if let Some(cr) = aux_string(&record, TAG_CR) {
            if cr != barcode_seq(&cb) {
                m.barcode_corrected_reads += 1;
            }
        }

        // UMI correction: encode(UR) vs encode(UB).
        if let Some(ur) = aux_string(&record, TAG_UR) {
            if let Ok(ur_enc) = barcode_str_to_u64(&ur) {
                if ur_enc != umi {
                    m.umi_corrected_reads += 1;
                }
            }
        }

        // Mapping classification.
        if record.tid() < 0 {
            m.unmapped_reads += 1;
        } else if record.mapq() != CONF_MAPPED_MAPQ {
            m.nonconf_mapped_reads += 1;
        } else if has_aux(&record, TAG_TX) {
            // tid != -1 && mapq == 255 && transcriptome-compatible
            m.conf_mapped += 1;
        }
        // else: uniquely mapped but no TX tag -> counts in `reads` only.
    }

    if verbose {
        eprintln!("Finished BAM parsing: {n} records total.");
    }
    Ok(metrics)
}

// ---------------------------------------------------------------------------
// File location helpers
// ---------------------------------------------------------------------------

/// Find a file by name in `dir`, then in `dir/outs`.
fn locate(dir: &Path, name: &str) -> Option<PathBuf> {
    let direct = dir.join(name);
    if direct.exists() {
        return Some(direct);
    }
    let outs = dir.join("outs").join(name);
    if outs.exists() {
        return Some(outs);
    }
    None
}

// ---------------------------------------------------------------------------
// main — implementation of `mkmolinfo::main`
// ---------------------------------------------------------------------------

fn main() -> Result<()> {
    let cli = Cli::parse();

    let bam_path = locate(&cli.out_directory, BAM_FILENAME).ok_or_else(|| {
        anyhow!(
            "Could not find {BAM_FILENAME} in {}",
            cli.out_directory.display()
        )
    })?;
    let input_h5 = locate(&cli.out_directory, MOLECULE_INFO_FILENAME).ok_or_else(|| {
        anyhow!(
            "Could not find {MOLECULE_INFO_FILENAME} in {}",
            cli.out_directory.display()
        )
    })?;

    // 1. Copy the input h5 to the output path, then open the copy read-write.
    //    In --self-check mode we never write, so open the input file read-only and
    //    leave the output untouched.
    let file = if cli.self_check {
        hdf5::File::open(&input_h5).context("Could not open hdf5 file")?
    } else {
        std::fs::copy(&input_h5, &cli.output).context("Could not copy input h5 file")?;
        hdf5::File::open_rw(&cli.output)
            .context("Could not open hdf5 output file (copy of input).")?
    };

    // 2. Read the molecule list from the copied file (3 datasets, dataset order).
    let barcode_idx =
        read_u64_dataset(&file, "/barcode_idx", "h5 file did not contain barcode_idx")?;
    let umis = read_u64_dataset(&file, "/umi", "H5 file did not contain umi")?;
    let barcode_strings = read_barcode_strings(&file)?;

    if barcode_idx.len() != umis.len() {
        bail!(
            "barcode_idx ({}) and umi ({}) lengths differ",
            barcode_idx.len(),
            umis.len()
        );
    }
    let n_mol = barcode_idx.len();

    // 3. Resolve each molecule's barcode string -> 2-bit code; enforce a single
    //    GEM group (error path).
    let mut bc_codes: Vec<u64> = Vec::with_capacity(n_mol);
    let mut single_gem: Option<u16> = None;
    for &bi in &barcode_idx {
        let bc_str = barcode_strings
            .get(bi as usize)
            .ok_or_else(|| anyhow!("barcode_idx {bi} out of range"))?;
        let gg = gem_group_of(bc_str).unwrap_or(1);
        match single_gem {
            None => single_gem = Some(gg),
            Some(g) if g != gg => bail!(
                "Error: Barcode is from multiple GEM groups, only one is supported. \
                 Please email the author to enable this. BC: {bc_str}"
            ),
            _ => {}
        }
        bc_codes.push(barcode_str_to_u64(bc_str)?);
    }

    // 4. Build the set of valid (barcode, umi) keys.
    let valid_keys: HashSet<MolKey> = (0..n_mol).map(|i| (bc_codes[i], umis[i])).collect();

    // 5. Parse the BAM and accumulate per-molecule metrics.
    let metrics = parse_bam(&bam_path, &valid_keys, cli.verbose)?;

    // Sanity report (sanity report).
    let found = metrics.len();
    eprintln!(
        "Barcodes in molecule_info.h5 found in BAM = {found}, not found = {}",
        valid_keys.len().saturating_sub(found)
    );

    if cli.self_check {
        let total_reads: u64 = metrics.values().map(|m| m.reads as u64).sum();
        let total_conf: u64 = metrics.values().map(|m| m.conf_mapped as u64).sum();
        println!("self-check: molecules in file       = {n_mol}");
        println!("self-check: molecules joined in BAM  = {found}");
        println!("self-check: total joined reads       = {total_reads}");
        println!("self-check: total conf_mapped reads  = {total_conf}");
        if found == 0 {
            eprintln!(
                "self-check WARNING: 0 molecules joined — the file's /umi encoding may \
                 not match barcode_str_to_u64 (A=0,C=1,G=2,T=3, MSB-first), or CB/UB \
                 tags are absent from the BAM."
            );
        }
        return Ok(());
    }

    // 6. Assemble per-molecule output arrays in molecule order and write the
    //    seven output datasets: six u32 counts + the u64 barcode array.
    let mut reads = Vec::with_capacity(n_mol);
    let mut bc_corr = Vec::with_capacity(n_mol);
    let mut conf = Vec::with_capacity(n_mol);
    let mut nonconf = Vec::with_capacity(n_mol);
    let mut umi_corr = Vec::with_capacity(n_mol);
    let mut unmapped = Vec::with_capacity(n_mol);

    for i in 0..n_mol {
        let m = metrics
            .get(&(bc_codes[i], umis[i]))
            .cloned()
            .unwrap_or_default();
        reads.push(m.reads);
        bc_corr.push(m.barcode_corrected_reads);
        conf.push(m.conf_mapped);
        nonconf.push(m.nonconf_mapped_reads);
        umi_corr.push(m.umi_corrected_reads);
        unmapped.push(m.unmapped_reads);
    }

    // Dataset names defined for this utility. Note: "conf_mapped" (not "..._reads").
    write_array(&file, "reads", &Array1::from(reads))?;
    write_array(&file, "barcode_corrected_reads", &Array1::from(bc_corr))?;
    write_array(&file, "conf_mapped", &Array1::from(conf))?;
    write_array(&file, "nonconf_mapped_reads", &Array1::from(nonconf))?;
    write_array(&file, "umi_corrected_reads", &Array1::from(umi_corr))?;
    write_array(&file, "unmapped_reads", &Array1::from(unmapped))?;
    // 7th dataset: barcode as the 2-bit u64 encoding (this command writes this one
    // directly via DatasetBuilder::create + Container::write, separate from
    // write_array because the element type is u64 rather than u32).
    write_array(&file, "barcode", &Array1::from(bc_codes))?;

    eprintln!("Wrote molecule info to {}", cli.output.display());
    Ok(())
}
