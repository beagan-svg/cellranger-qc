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

## Release checklist

1. Update the package version in `pyproject.toml`.
2. Update `src/cellranger_qc/__init__.py`.
3. Add release notes to `CHANGELOG.md`.
4. Commit the release, create a `vX.Y.Z` tag, and push the tag.
