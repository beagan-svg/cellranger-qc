#!/usr/bin/env python3
"""Compare reconstructed molecule_info datasets against expected values."""

import json
import sys
from pathlib import Path

import h5py
import numpy as np

DATASETS = [
    "reads",
    "conf_mapped",
    "nonconf_mapped_reads",
    "unmapped_reads",
    "barcode_corrected_reads",
    "umi_corrected_reads",
    "barcode",
]


def flattened_or_none(h5_file: h5py.File, dataset: str) -> np.ndarray | None:
    """Return a flattened dataset array, or None if the dataset is absent."""
    if dataset not in h5_file:
        return None
    return np.asarray(h5_file[dataset][:]).ravel()


def compare(orig_path: Path, mine_path: Path, expected_path: Path) -> bool:
    """Print a dataset comparison table and return True when all datasets match."""
    expected = json.loads(expected_path.read_text())
    ok = True

    with h5py.File(orig_path, "r") as orig, h5py.File(mine_path, "r") as mine:
        print(f"{'dataset':<26} {'original':<22} {'reconstruction':<22} {'expected':<22} match")
        print("-" * 100)
        for dataset in DATASETS:
            original_values = flattened_or_none(orig, dataset)
            reconstructed_values = flattened_or_none(mine, dataset)
            expected_values = np.array(expected[dataset])

            same_original_reconstructed = (
                original_values is not None
                and reconstructed_values is not None
                and np.array_equal(original_values, reconstructed_values)
            )
            same_original_expected = original_values is not None and np.array_equal(
                original_values, expected_values
            )
            same_reconstructed_expected = reconstructed_values is not None and np.array_equal(
                reconstructed_values, expected_values
            )
            row_ok = (
                same_original_reconstructed
                and same_original_expected
                and same_reconstructed_expected
            )
            ok = ok and row_ok
            match_label = (
                "YES"
                if row_ok
                else (
                    "NO  "
                    f"(orig==mine:{same_original_reconstructed} "
                    f"orig==exp:{same_original_expected})"
                )
            )
            original_list = None if original_values is None else list(original_values)
            reconstructed_list = (
                None if reconstructed_values is None else list(reconstructed_values)
            )
            print(
                f"{dataset:<26} {str(original_list):<22} "
                f"{str(reconstructed_list):<22} {str(list(expected_values)):<22} "
                f"{match_label}"
            )

    print("-" * 100)
    print("ALL DATASETS IDENTICAL (original == reconstruction == expected):", ok)
    return ok


def main() -> int:
    """CLI entry point."""
    if len(sys.argv) != 4:
        print("usage: compare_h5.py ORIGINAL_H5 RECONSTRUCTED_H5 EXPECTED_JSON", file=sys.stderr)
        return 2
    return 0 if compare(*(Path(value) for value in sys.argv[1:])) else 1


if __name__ == "__main__":
    raise SystemExit(main())
