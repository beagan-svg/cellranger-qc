# cellranger-qc

Post-run quality-control utilities for Cell Ranger GEX, multiome, and ATAC outputs.

`cellranger-qc` packages two production command-line workflows:

- `cellranger-gex-qc` computes post-alignment GEX/multiome QC tables, doublet scores, intron/exon matrices, count matrices, and library-level summaries.
- `cellranger-atac-qc` computes ATAC TSS enrichment, promoter ratio, nucleosome metrics, merges the metrics into a GEX AnnData object, and exports a CSV.

## Requirements

- Python 3.10, 3.11, or 3.12
- [`uv`](https://docs.astral.sh/uv/) for environment and dependency management
- Cell Ranger output directories with the expected `outs/` files
- For ATAC QC, a bgzip-compressed and tabix-indexed fragments file readable by `pysam.TabixFile`

Python 3.13+ is not enabled because several scientific Python and bioinformatics dependencies can lag new Python releases.

## Installation

Install the runtime environment:

```bash
uv sync
```

Install the development environment:

```bash
uv sync --extra dev
```

Confirm the CLIs are available:

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

`--libs` must point to a CSV file with one row per library.

| Column | Description |
| --- | --- |
| `ar_dir` | Parent directory containing the alignment directory |
| `ar_id` | Alignment ID; appended to `ar_dir` to find the Cell Ranger run |
| `library_prep` | Library prep name used in sample IDs and output names |
| `cell_prep_type` | `Cells` uses a 1,500 gene threshold; other values use 1,000 |
| `load_name` | Prefix used to create `cell_member` |
| `alignment_method` | Used to identify Cell Ranger ARC/multi outputs |
| `expc_cell_capture` | Expected cell capture count for usable-cell percentage |

### Expected Inputs

For each library, the command expects these Cell Ranger files:

- `outs/filtered_feature_bc_matrix/matrix.mtx.gz`
- `outs/filtered_feature_bc_matrix/features.tsv.gz`
- `outs/filtered_feature_bc_matrix/barcodes.tsv.gz`
- `outs/*molecule_info.h5`
- `outs/*summary.csv`
- `outs/web_summary.html`
- Either `outs/per_barcode_metrics.csv` or `outs/*molecule_info_new.h5`

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
  --load-name LOAD_001 \
  --library-prep-name LIB_PREP_001 \
  --alignment-fs-id AR123 \
  --annotation-file genes.gtf.gz \
  --atac-fragments-path fragments.tsv.gz \
  --per-barcode-metrics-path per_barcode_metrics.csv
```

### Expected Inputs

`--output-path` should contain:

- `gex_<alignment_fs_id>.h5ad`
- `samp.dat_<alignment_fs_id>.csv`

The `samp.dat` file may use either `bc` or `barcodes` for the raw barcode column.

Additional required inputs:

- `--annotation-file`: GTF annotation file, optionally gzipped
- `--atac-fragments-path`: bgzip-compressed, tabix-indexed ATAC fragments file
- `--per-barcode-metrics-path`: Cell Ranger `per_barcode_metrics.csv`

### Outputs

The ATAC workflow:

- Updates `gex_<alignment_fs_id>.h5ad` by joining ATAC metrics into `obs`
- Writes `atac_qc_<alignment_fs_id>.csv`

## Operational Notes

- Commands log progress to standard output with timestamps and step summaries.
- GEX doublet scoring uses a deterministic random seed for synthetic doublet generation.
- ATAC insertion counting parallelizes chromosome batches internally.
- Large fragment files should be stored on fast local or high-throughput shared storage.
- Generated matrices, HDF5 files, CSVs, and PDFs are ignored by git by default.

## Development

Install development dependencies:

```bash
uv sync --extra dev
```

Run checks:

```bash
uv run ruff check .
uv run pytest
```

GitHub Actions runs linting and tests on Python 3.10, 3.11, and 3.12.

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
