# cellranger-qc

Post-run quality-control utilities for Cell Ranger outputs.

This package currently provides two command-line workflows:

- `cellranger-gex-qc`: post-alignment GEX/multiome QC, including per-cell gene counts, doublet scores, intron/exon matrices, count matrices, and OCS summary output.
- `cellranger-atac-qc`: ATAC QC for Cell Ranger ARC/multiome outputs, including TSS enrichment, promoter ratio, nucleosome metrics, AnnData metadata merge, and CSV export.

## Requirements

- Python 3.10, 3.11, or 3.12
- Cell Ranger output directories with the expected `outs/` files
- For ATAC QC, a bgzip-compressed and tabix-indexed fragments file readable by `pysam.TabixFile`

Python 3.13+ is intentionally not enabled yet because several scientific Python and bioinformatics dependencies can lag new Python releases.

## Install

Clone and install in editable mode:

```bash
git clone git@github.com:beagan-svg/cellranger-qc.git
cd cellranger-qc
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

For runtime-only installation:

```bash
python -m pip install .
```

Check the installed commands:

```bash
cellranger-gex-qc --help
cellranger-atac-qc --help
```

## GEX / Multiome QC

Run:

```bash
cellranger-gex-qc \
  --libs libs.csv \
  --out-dir qc_output \
  --num-cores 16
```

The `libs.csv` file is expected to contain one row per library. The script uses these columns:

| Column | Description |
| --- | --- |
| `ar_dir` | Parent directory containing the alignment directory |
| `ar_id` | Alignment ID; appended to `ar_dir` to find the Cell Ranger run |
| `library_prep` | Library prep name used in sample IDs and output names |
| `cell_prep_type` | `Cells` uses a 1,500 gene threshold; other values use 1,000 |
| `load_name` | Prefix used to create `cell_member` |
| `alignment_method` | Used to identify Cell Ranger ARC/multi outputs |
| `expc_cell_capture` | Expected cell capture count for usable-cell percentage |

Expected Cell Ranger files include:

- `outs/filtered_feature_bc_matrix/matrix.mtx.gz`
- `outs/filtered_feature_bc_matrix/features.tsv.gz`
- `outs/filtered_feature_bc_matrix/barcodes.tsv.gz`
- `outs/*molecule_info.h5`
- `outs/*summary.csv`
- `outs/web_summary.html`
- Either `outs/per_barcode_metrics.csv` or `outs/*molecule_info_new.h5`

Main outputs:

- `matrix/count_<ar_id>.mtx`
- `matrix/intron_<ar_id>.mtx`
- `matrix/exon_<ar_id>.mtx`
- `samp_dat_<ar_id>.csv`
- `<library_prep>.doubscore.pdf`
- `ocs_summary.csv`

## ATAC QC

Run:

```bash
cellranger-atac-qc \
  --output-path qc_output \
  --load-name LOAD_001 \
  --library-prep-name LIB_PREP_001 \
  --alignment-fs-id AR123 \
  --annotation-file genes.gtf.gz \
  --atac-fragments-path fragments.tsv.gz \
  --per-barcode-metrics-path per_barcode_metrics.csv
```

The output directory should contain:

- `gex_<alignment_fs_id>.h5ad`
- `samp.dat_<alignment_fs_id>.csv`

The `samp.dat` file may use either `bc` or `barcodes` for the raw barcode column.

Main outputs:

- Updates `gex_<alignment_fs_id>.h5ad` by joining ATAC metrics into `obs`
- Writes `atac_qc_<alignment_fs_id>.csv`

## Development

Install development dependencies:

```bash
python -m pip install -e ".[dev]"
```

Run checks:

```bash
ruff check .
pytest
```

The GitHub Actions workflow runs linting and tests on Python 3.10, 3.11, and 3.12.

## Versioning

This project uses semantic versioning:

- `MAJOR`: incompatible CLI, input, or output changes
- `MINOR`: new features that remain backward compatible
- `PATCH`: bug fixes and small documentation updates

When releasing:

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

## Repository Setup

This repository is intended to publish to:

```bash
git@github.com:beagan-svg/cellranger-qc.git
```

If starting from a fresh local checkout:

```bash
git init
git remote add origin git@github.com:beagan-svg/cellranger-qc.git
git branch -M main
git add .
git commit -m "Initial production package"
git push -u origin main
```
