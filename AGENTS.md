# Agent Instructions

These instructions apply to the whole repository.

## Project Shape

- Python package source lives in `src/cellranger_qc`.
- Python tests live in `tests`.
- The Rust molecule-info helper lives in `mkmolinfo-rs`.
- Generated Cell Ranger outputs, HDF5 files, matrices, PDFs, and Rust build
  artifacts should not be committed.

## Checks

Run Python checks before handing off changes:

```bash
uv run ruff check .
uv run pytest
```

Run Rust checks when touching `mkmolinfo-rs`:

```bash
cargo fmt --check --manifest-path mkmolinfo-rs/Cargo.toml
cargo test --manifest-path mkmolinfo-rs/Cargo.toml
```

## Style

- Prefer clear scientific naming over terse abbreviations.
- Keep public CLI behavior and output column names stable unless the changelog
  calls out a breaking change.
- Document biological assumptions in docstrings when they affect thresholds,
  bins, or QC interpretation.
- Keep generated data out of tests; use small synthetic fixtures or mocks.
