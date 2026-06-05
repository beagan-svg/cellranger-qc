# mkmolinfo parity test

Verifies the Rust reconstruction produces output identical to the original
`mkmolinfo` binary. The original is a Linux x86-64 build, so it runs under
`docker --platform linux/amd64` (statically linked — only needs glibc); the
reconstruction runs natively on the host.

## Run

```bash
./run_parity.sh
```

Prints a table comparing all 7 datasets and ends with
`ALL DATASETS IDENTICAL (original == reconstruction == expected): True`.

## Requirements

- Docker (with amd64 emulation)
- HDF5 1.10 or 1.14 and `htslib` (to build/run the reconstruction). The
  current Rust HDF5 binding does not accept Homebrew HDF5 2.x.
- The original `mkmolinfo` binary in the project root

## Files

- `gen_testdata.py` — writes `testdata/outs/{molecule_info.h5, possorted_genome_bam.bam}`
  plus `expected.json`. The reads exercise every counter (conf/nonconf/unmapped,
  barcode/UMI correction), the mapq-255-without-`TX` edge case, a multi-read
  molecule, a zero-read molecule, and a decoy not in the molecule set.
- `compare_h5.py` — diffs the 7 datasets of two outputs against `expected.json`.
