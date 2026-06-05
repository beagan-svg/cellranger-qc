#!/usr/bin/env bash
# Validation test: run the Rust implementation on synthetic data, then compare
# all 7 output datasets against expected values. See README.md.
set -euo pipefail

test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rs_dir="$(dirname "$test_dir")"
root_dir="$(dirname "$rs_dir")"

hdf5_dir="$(brew --prefix hdf5@1.10)"
htslib_dir="$(brew --prefix htslib)"
export HDF5_DIR="$hdf5_dir" HTSLIB_DIR="$htslib_dir"
export PKG_CONFIG_PATH="$htslib_dir/lib/pkgconfig:$hdf5_dir/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export RUSTFLAGS="-C link-args=-Wl,-rpath,$hdf5_dir/lib"
export DYLD_FALLBACK_LIBRARY_PATH="$hdf5_dir/lib:$htslib_dir/lib"

echo "==> Building implementation (host)"
cargo build --quiet --manifest-path "$rs_dir/Cargo.toml"

echo "==> Generating synthetic data"
docker run --rm --platform linux/amd64 -v "$root_dir":/work -w /work/mkmolinfo-rs/test python:3.11-slim bash -c '
  set -e
  pip install --quiet --no-cache-dir h5py numpy pysam
  rm -rf testdata && mkdir -p testdata/outs
  python gen_testdata.py testdata/outs >/dev/null
'

echo "==> Running IMPLEMENTATION (host)"
"$rs_dir/target/debug/mkmolinfo" "$test_dir/testdata" -o "$test_dir/out_mine.h5"

echo "==> Comparing all 7 datasets (docker linux/amd64)"
docker run --rm --platform linux/amd64 -v "$root_dir":/work -w /work/mkmolinfo-rs/test python:3.11-slim bash -c '
  pip install --quiet --no-cache-dir h5py numpy
  python compare_h5.py out_mine.h5 testdata/outs/expected.json
'
