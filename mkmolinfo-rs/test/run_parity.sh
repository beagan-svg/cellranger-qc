#!/usr/bin/env bash
# Parity test: run the original mkmolinfo binary (Linux, via docker linux/amd64)
# and the Rust reconstruction (host) on the same synthetic data, then diff all 7
# output datasets. See README.md.
set -euo pipefail

test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rs_dir="$(dirname "$test_dir")"
root_dir="$(dirname "$rs_dir")"   # holds the original `mkmolinfo` binary

hdf5_dir="$(brew --prefix hdf5@1.10)"
htslib_dir="$(brew --prefix htslib)"
export HDF5_DIR="$hdf5_dir" HTSLIB_DIR="$htslib_dir"
export PKG_CONFIG_PATH="$htslib_dir/lib/pkgconfig:$hdf5_dir/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export RUSTFLAGS="-C link-args=-Wl,-rpath,$hdf5_dir/lib"
export DYLD_FALLBACK_LIBRARY_PATH="$hdf5_dir/lib:$htslib_dir/lib"

echo "==> Building reconstruction (host)"
cargo build --quiet --manifest-path "$rs_dir/Cargo.toml"

echo "==> Generating data + running ORIGINAL binary (docker linux/amd64)"
docker run --rm --platform linux/amd64 -v "$root_dir":/work -w /work/mkmolinfo-rs/test python:3.11-slim bash -c '
  set -e
  pip install --quiet --no-cache-dir h5py numpy pysam
  chmod +x /work/mkmolinfo
  rm -rf testdata && mkdir -p testdata/outs
  python gen_testdata.py testdata/outs >/dev/null
  /work/mkmolinfo testdata -o out_orig.h5
'

echo "==> Running RECONSTRUCTION (host)"
"$rs_dir/target/debug/mkmolinfo" "$test_dir/testdata" -o "$test_dir/out_mine.h5"

echo "==> Comparing all 7 datasets (docker linux/amd64)"
docker run --rm --platform linux/amd64 -v "$root_dir":/work -w /work/mkmolinfo-rs/test python:3.11-slim bash -c '
  pip install --quiet --no-cache-dir h5py numpy
  python compare_h5.py out_orig.h5 out_mine.h5 testdata/outs/expected.json
'
