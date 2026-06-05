# molecule_info_new.h5 helper validation test

Verifies the helper against a small synthetic dataset. Docker is used to run
the comparison in a consistent Linux environment; the helper runs natively on
the host.

## Run

```bash
./run_validation.sh
```

Prints a table comparing all 7 datasets and ends with
`ALL DATASETS IDENTICAL (generated == expected): True`.

## Requirements

- Docker (with amd64 emulation)
- HDF5 1.10 or 1.14 and `htslib` (to build/run the helper). The
  current Rust HDF5 binding does not accept Homebrew HDF5 2.x.

## Files

- `gen_testdata.py` — writes `testdata/outs/{molecule_info.h5, possorted_genome_bam.bam}`
  plus `expected.json`. The reads exercise every counter (conf/nonconf/unmapped,
  barcode/UMI correction), the mapq-255-without-`TX` edge case, a multi-read
  molecule, a zero-read molecule, and a decoy not in the molecule set.
- `compare_h5.py` — compares the generated datasets against `expected.json`.
