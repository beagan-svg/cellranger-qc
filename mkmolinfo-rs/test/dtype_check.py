#!/usr/bin/env python3
"""Compare generated molecule_info dataset dtypes."""

from pathlib import Path

import h5py

DATASETS = [
    "reads",
    "conf_mapped",
    "nonconf_mapped_reads",
    "unmapped_reads",
    "barcode_corrected_reads",
    "umi_corrected_reads",
    "barcode",
]
COPIED_THROUGH_DATASETS = ["barcode_idx", "umi", "barcodes"]


def main() -> int:
    """CLI entry point."""
    generated_path = Path("out_mine.h5")

    with h5py.File(generated_path) as generated:
        print(f"{'dataset':<26}{'dtype':<14}expected")
        all_ok = True
        for dataset in DATASETS:
            generated_dtype = generated[dataset].dtype
            expected_dtype = "uint64" if dataset == "barcode" else "uint32"
            dtype_matches = str(generated_dtype) == expected_dtype
            all_ok = all_ok and dtype_matches
            print(f"{dataset:<26}{str(generated_dtype):<14}{expected_dtype}")
        copied_through_preserved = all(dataset in generated for dataset in COPIED_THROUGH_DATASETS)
        print("DTYPES IDENTICAL:", all_ok)
        print("copied-through preserved:", copied_through_preserved)
    return 0 if all_ok and copied_through_preserved else 1


if __name__ == "__main__":
    raise SystemExit(main())
