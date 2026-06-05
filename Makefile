.PHONY: sync lint test check rust-fmt rust-test rust-check

sync:
	uv sync --extra dev

lint:
	uv run ruff check .

test:
	uv run pytest

check: lint test

rust-fmt:
	cargo fmt --check --manifest-path mkmolinfo-rs/Cargo.toml

rust-test:
	cargo test --manifest-path mkmolinfo-rs/Cargo.toml

rust-check: rust-fmt rust-test
