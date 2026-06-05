# Contributing

## Local setup

```bash
uv sync --extra dev
```

## Checks

Run these before opening a pull request:

```bash
uv run ruff check .
uv run pytest
```

If you changed `mkmolinfo-rs`, also run:

```bash
cargo fmt --check --manifest-path mkmolinfo-rs/Cargo.toml
cargo test --manifest-path mkmolinfo-rs/Cargo.toml
```

The same commands are available as Makefile shortcuts:

```bash
make check
make rust-check
```

## Release checklist

1. Update the package version in `pyproject.toml`.
2. Update `src/cellranger_qc/__init__.py`.
3. Add release notes to `CHANGELOG.md`.
4. Commit the release, create a `vX.Y.Z` tag, and push the tag.
