# mkmolinfo-rs

`mkmolinfo-rs` is a Rust reconstruction of the legacy 10x Genomics
`mkmolinfo` utility. It converts a Cell Ranger output directory containing a
new-style `molecule_info.h5` and `possorted_genome_bam.bam` into the older
`molecule_info_new.h5` shape expected by parts of the GEX QC workflow.

## Build

The HDF5 Rust binding used by this crate does not currently accept Homebrew
HDF5 2.x. Use HDF5 1.10 or 1.14 when building locally on macOS.

```bash
cargo build --manifest-path mkmolinfo-rs/Cargo.toml --release
```

The release binary is written to:

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

The parity test harness under `mkmolinfo-rs/test` compares the reconstruction
against the original binary when Docker and the original `mkmolinfo` executable
are available.
