# cellranger-qc

[![CI](https://github.com/beagan-svg/cellranger-qc/actions/workflows/ci.yml/badge.svg)](https://github.com/beagan-svg/cellranger-qc/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.12-blue)](https://www.python.org/)

Command-line tools for calculating quality-control metrics from Cell Ranger GEX,
multiome, and ATAC outputs.

The package provides two commands:

- `cellranger-gex-qc` creates per-cell GEX metrics, count matrices, doublet scores,
  and a library summary.
- `cellranger-atac-qc` calculates TSS enrichment, promoter coverage, fragment counts,
  and nucleosome metrics.

## Requirements

- Python 3.12
- Cell Ranger output files for the workflow being processed
- A bgzip-compressed, tabix-indexed fragments file for ATAC QC
- Rust stable only when building the optional `mkmolinfo` helper

## Installation

Clone the repository and install the locked dependencies:

```bash
git clone https://github.com/beagan-svg/cellranger-qc.git
cd cellranger-qc
uv sync --frozen
```

Confirm that both commands are available:

```bash
uv run cellranger-gex-qc --help
uv run cellranger-atac-qc --help
```

## GEX and multiome QC

Create a CSV manifest with one row per library:

```csv
cellranger_run_dir,library_prep,cell_prep_type,load_name,alignment_method,expc_cell_capture
/data/run-001,run-001,Cells,load-001,CELL_RANGER_COUNT,5000
```

The manifest fields are:

| Field | Description |
| --- | --- |
| `cellranger_run_dir` | Cell Ranger run directory containing `outs/` |
| `library_prep` | Name included in output filenames and cell IDs |
| `cell_prep_type` | `Cells` uses the 1,500-gene filter; other values use 1,000 |
| `load_name` | Prefix used for the `cell_member` value |
| `alignment_method` | Cell Ranger workflow identifier |
| `expc_cell_capture` | Expected cell count used to calculate usable yield |

Run the workflow:

```bash
uv run cellranger-gex-qc \
  --libs libraries.csv \
  --out-dir qc-output \
  --num-cores 16
```

Each Cell Ranger `outs/` directory must contain:

- `filtered_feature_bc_matrix/matrix.mtx.gz`
- `filtered_feature_bc_matrix/features.tsv.gz`
- `filtered_feature_bc_matrix/barcodes.tsv.gz`
- `molecule_info.h5`
- a Cell Ranger summary CSV
- `web_summary.html`
- either `per_barcode_metrics.csv` or `molecule_info_new.h5`

The command writes:

```text
qc-output/
├── <library_prep>.doubscore.pdf
├── matrix/
│   ├── count_<library_prep>.mtx
│   ├── exon_<library_prep>.mtx
│   └── intron_<library_prep>.mtx
├── ocs_summary.csv
└── samp_dat_<library_prep>.csv
```

For Cell Ranger count runs that do not already provide `molecule_info_new.h5`,
build and run the optional Rust helper described in
[`docs/mkmolinfo.md`](docs/mkmolinfo.md).

## ATAC QC

ATAC QC needs a GTF annotation, `per_barcode_metrics.csv`, and a fragments file
with its sibling `.tbi` index.

```bash
uv run cellranger-atac-qc \
  --output-path qc-output \
  --annotation-file /data/reference/genes.gtf.gz \
  --atac-fragments-path /data/run-001/outs/atac_fragments.tsv.gz \
  --per-barcode-metrics-path /data/run-001/outs/per_barcode_metrics.csv
```

The command writes `qc-output/atac_qc.csv`. The file contains passing cell
barcodes and their TSS enrichment, promoter ratio, fragment counts, and
nucleosome measurements.

## Development

Development commands are defined in the `justfile`. Install the environment and
run the Python checks with:

```bash
just sync
just check
```

When changing the Rust helper, run:

```bash
just rust-check
```
