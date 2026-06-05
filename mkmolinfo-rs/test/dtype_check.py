#!/usr/bin/env python3
"""Compare reconstructed molecule_info dataset dtypes."""

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
    original_path = Path("out_orig.h5")
    reconstructed_path = Path("out_mine.h5")

    with h5py.File(original_path) as original, h5py.File(reconstructed_path) as reconstructed:
        print(f"{'dataset':<26}{'orig dtype':<14}{'mine dtype':<14}match")
        all_ok = True
        for dataset in DATASETS:
            original_dtype = original[dataset].dtype
            reconstructed_dtype = reconstructed[dataset].dtype
            dtype_matches = original_dtype == reconstructed_dtype
            all_ok = all_ok and dtype_matches
            print(
                f"{dataset:<26}{str(original_dtype):<14}"
                f"{str(reconstructed_dtype):<14}{dtype_matches}"
            )
        copied_through_preserved = all(
            dataset in original and dataset in reconstructed for dataset in COPIED_THROUGH_DATASETS
        )
        print("DTYPES IDENTICAL:", all_ok)
        print("copied-through preserved:", copied_through_preserved)
    return 0 if all_ok and copied_through_preserved else 1


if __name__ == "__main__":
    raise SystemExit(main())
