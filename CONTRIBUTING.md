# Contributing

## Local setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

## Checks

Run these before opening a pull request:

```bash
ruff check .
pytest
```

## Release checklist

1. Update the package version in `pyproject.toml`.
2. Update `src/cellranger_qc/__init__.py`.
3. Add release notes to `CHANGELOG.md`.
4. Commit the release, create a `vX.Y.Z` tag, and push the tag.
