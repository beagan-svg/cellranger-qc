# mkmolinfo-rs

`mkmolinfo-rs` is a Rust implementation of the 10x Genomics
`mkmolinfo` utility. It reads a Cell Ranger output directory containing
`molecule_info.h5` and `possorted_genome_bam.bam`, then writes
`molecule_info_new.h5` for GEX QC workflows that need that file.

## Build

The HDF5 Rust binding used by this crate does not currently accept Homebrew
HDF5 2.x. Use HDF5 1.10 or 1.14 when building locally on macOS.

```bash
cargo build --manifest-path mkmolinfo-rs/Cargo.toml --release
```

The release executable is written to:

```text
mkmolinfo-rs/target/release/mkmolinfo
```

## Run

```bash
mkmolinfo-rs/target/release/mkmolinfo /path/to/cellranger/outs
```

By default the command writes `molecule_info_new.h5` in the current working
directory. Use `--output` to choose a different path:

```bash
mkmolinfo-rs/target/release/mkmolinfo /path/to/cellranger/outs \
  --output /path/to/molecule_info_new.h5
```

## Checks

```bash
cargo fmt --check --manifest-path mkmolinfo-rs/Cargo.toml
cargo test --manifest-path mkmolinfo-rs/Cargo.toml
```

The validation harness under `mkmolinfo-rs/test` compares generated datasets
against expected values when Docker is available.
