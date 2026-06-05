# cellranger-qc

[![CI](https://github.com/beagan-svg/cellranger-qc/actions/workflows/ci.yml/badge.svg)](https://github.com/beagan-svg/cellranger-qc/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org/)
[![uv](https://img.shields.io/badge/package%20manager-uv-6340ac)](https://docs.astral.sh/uv/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Post-run quality-control tools for Cell Ranger GEX, multiome, and ATAC outputs.

This repository provides:

- `cellranger-gex-qc`: GEX/multiome post-alignment QC, doublet scores, intron/exon matrices, count matrices, and library summaries.
- `cellranger-atac-qc`: portable ATAC QC from fragments, per-barcode metrics, and a GTF annotation.
- `mkmolinfo-rs`: optional Rust helper for reconstructing legacy `molecule_info_new.h5` files used by older GEX QC workflows.

## Requirements

- Python 3.10, 3.11, or 3.12
- [`uv`](https://docs.astral.sh/uv/) for Python environment management
- Cell Ranger output files for the workflow you are running
- Rust stable only if you build `mkmolinfo-rs`

Python 3.13+ is not enabled yet because some scientific Python and bioinformatics dependencies can lag new Python releases.

## Setup

Install runtime dependencies:

```bash
uv sync
```

Install development dependencies:

```bash
uv sync --extra dev
```

Check the installed commands:

```bash
uv run cellranger-gex-qc --help
uv run cellranger-atac-qc --help
```

## GEX / Multiome QC

Run:

```bash
uv run cellranger-gex-qc \
  --libs libs.csv \
  --out-dir qc_output \
  --num-cores 16
```

### Library Manifest

`--libs` is a CSV with one row per library.

| Column | Description |
| --- | --- |
| `ar_dir` | Parent directory containing the alignment directory |
| `ar_id` | Alignment ID appended to `ar_dir` to find the Cell Ranger run |
| `library_prep` | Library prep name used in sample IDs and output names |
| `cell_prep_type` | `Cells` uses a 1,500 gene threshold; other values use 1,000 |
| `load_name` | Prefix used to create `cell_member` |
| `alignment_method` | Used to identify Cell Ranger ARC/multi outputs |
| `expc_cell_capture` | Expected cell capture count for usable-cell percentage |

### Expected Inputs

For each manifest row, the workflow expects these files under the resolved Cell Ranger `outs` directory:

- `filtered_feature_bc_matrix/matrix.mtx.gz`
- `filtered_feature_bc_matrix/features.tsv.gz`
- `filtered_feature_bc_matrix/barcodes.tsv.gz`
- `*molecule_info.h5`
- `*summary.csv`
- `web_summary.html`
- Either `per_barcode_metrics.csv` or `*molecule_info_new.h5`

### Outputs

The GEX/multiome workflow writes:

- `matrix/count_<ar_id>.mtx`
- `matrix/intron_<ar_id>.mtx`
- `matrix/exon_<ar_id>.mtx`
- `samp_dat_<ar_id>.csv`
- `<library_prep>.doubscore.pdf`
- `ocs_summary.csv`

## ATAC QC

Run:

```bash
uv run cellranger-atac-qc \
  --output-path qc_output \
  --annotation-file genes.gtf.gz \
  --atac-fragments-path fragments.tsv.gz \
  --per-barcode-metrics-path per_barcode_metrics.csv
```

### Expected Inputs

- `--annotation-file`: GTF annotation file, optionally gzipped
- `--atac-fragments-path`: bgzip-compressed, tabix-indexed Cell Ranger ATAC `fragments.tsv.gz`
- `--per-barcode-metrics-path`: Cell Ranger `per_barcode_metrics.csv`
- `--output-path`: directory where `atac_qc.csv` will be written

The fragments file must have a sibling tabix index, usually `fragments.tsv.gz.tbi`.

### Outputs

The ATAC workflow writes `atac_qc.csv` with passing cells and `atac_`-prefixed metrics, including:

- `atac_tss_enrichment`
- `atac_reads_in_tss`
- `atac_n_frags`
- `atac_n_nucleosome_free_frags`
- `atac_n_mono_frags`
- `atac_n_di_frags`
- `atac_n_multi_frags`
- `atac_nucleosome_ratio`
- `atac_reads_in_promoter`
- `atac_promoter_ratio`

## Optional Legacy molecule_info Helper

`mkmolinfo-rs` reconstructs a legacy `molecule_info_new.h5` from a Cell Ranger output directory. Build it only if your GEX workflow needs that compatibility file.

```bash
cargo build --manifest-path mkmolinfo-rs/Cargo.toml --release
```

See [docs/mkmolinfo.md](docs/mkmolinfo.md) for usage, HDF5 version notes, and parity-test details.

## Operational Notes

- Commands log progress to standard output with timestamps and step summaries.
- GEX doublet scoring uses a deterministic random seed.
- ATAC insertion counting parallelizes chromosome batches internally.
- Large fragment files should live on fast local or high-throughput shared storage.
- Generated matrices, HDF5 files, CSVs, PDFs, Rust build outputs, and local environments are ignored by git.

## Development

Run Python checks:

```bash
uv sync --extra dev
uv run ruff check .
uv run pytest
```

Or use Makefile shortcuts:

```bash
make check
```

Run Rust checks if `mkmolinfo-rs` changed:

```bash
make rust-check
```

GitHub Actions runs Python linting/tests on Python 3.10, 3.11, and 3.12, plus Rust formatting/tests for `mkmolinfo-rs`.

## Versioning

This project uses semantic versioning:

- `MAJOR`: incompatible CLI, input, or output changes
- `MINOR`: backward-compatible features
- `PATCH`: backward-compatible fixes and documentation updates

Release checklist:

1. Update `version` in `pyproject.toml`.
2. Update `__version__` in `src/cellranger_qc/__init__.py`.
3. Add a new section to `CHANGELOG.md`.
4. Commit and tag the release:

```bash
git add pyproject.toml src/cellranger_qc/__init__.py CHANGELOG.md
git commit -m "Release v0.1.0"
git tag v0.1.0
git push origin main --tags
```

## License

This project is released under the [MIT License](LICENSE).
