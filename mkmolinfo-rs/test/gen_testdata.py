#!/usr/bin/env python3
"""Generate a small, fully controlled Cell Ranger-style input for validation testing.

Writes  <outdir>/molecule_info.h5   with /barcodes, /barcode_idx, /umi
and     <outdir>/possorted_genome_bam.bam  with reads carrying CB/UB/CR/UR/TX
tags and chosen mapq/unmapped state, so every per-molecule counter is exercised.

The /umi values are computed with the same 2-bit scheme used by
barcode_str_to_u64 (A=0,C=1,G=2,T=3, MSB-first).
"""

import sys
from pathlib import Path

import h5py
import numpy as np
import pysam


def enc(seq: str) -> int:
    m = {"A": 0, "C": 1, "G": 2, "T": 3}
    acc = 0
    for b in seq:
        acc = (acc << 2) | m[b]
    return acc


import json


def main() -> int:
    """Generate synthetic molecule_info and BAM fixtures."""
    if len(sys.argv) != 2:
        print("usage: gen_testdata.py OUT_DIRECTORY", file=sys.stderr)
        return 2

    outdir = Path(sys.argv[1])
    outdir.mkdir(parents=True, exist_ok=True)

    # Barcodes are 16 bp; UMIs are 10 bp.
    barcode_0 = "AAACCCGGGTTTAAAC"
    barcode_1 = "ACGTACGTACGTACGT"
    umi_0 = "AAAAAAAAAA"
    umi_1 = "ACACACACAC"
    umi_2 = "GGGGTTTTAA"

    barcodes = [barcode_0, barcode_1]  # /barcodes stores raw 16 bp strings, no GEM suffix.
    molecules = [
        (0, umi_0),  # M0 barcode_0, umi_0 -> many reads
        (0, umi_1),  # M1 barcode_0, umi_1 -> 1 read
        (1, umi_2),  # M2 barcode_1, umi_2 -> 2 reads
        (1, umi_0),  # M3 barcode_1, umi_0 -> no reads
    ]
    barcode_idx = np.array([molecule[0] for molecule in molecules], dtype=np.uint64)
    umi = np.array([enc(molecule[1]) for molecule in molecules], dtype=np.uint64)

    with h5py.File(outdir / "molecule_info.h5", "w") as h5_file:
        max_barcode_len = max(len(barcode) for barcode in barcodes)
        barcode_array = np.array(
            [barcode.encode("ascii") for barcode in barcodes],
            dtype=f"S{max_barcode_len}",
        )
        h5_file.create_dataset("barcodes", data=barcode_array)
        h5_file.create_dataset("barcode_idx", data=barcode_idx)
        h5_file.create_dataset("umi", data=umi)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1_000_000}],
    }

    # (qname, cb, ub, cr, ur, mapq, unmapped, has_tx)
    reads = [
        ("r1", f"{barcode_0}-1", umi_0, barcode_0, umi_0, 255, False, True),
        ("r2", f"{barcode_0}-1", umi_0, "AAACCCGGGTTTAAAG", umi_0, 255, False, True),
        ("r3", f"{barcode_0}-1", umi_0, barcode_0, "AAAAAAAAAT", 30, False, True),
        ("r4", f"{barcode_0}-1", umi_0, barcode_0, umi_0, 0, True, False),
        ("r5", f"{barcode_0}-1", umi_0, barcode_0, umi_0, 255, False, False),
        ("r6", f"{barcode_0}-1", umi_1, barcode_0, umi_1, 255, False, True),
        ("r7", f"{barcode_1}-1", umi_2, barcode_1, umi_2, 255, False, True),
        ("r8", f"{barcode_1}-1", umi_2, barcode_1, umi_2, 0, False, True),
        (
            "r9",
            "GGGGGGGGGGGGGGGG-1",
            "TTTTTTTTTT",
            "GGGGGGGGGGGGGGGG",
            "TTTTTTTTTT",
            255,
            False,
            True,
        ),
    ]

    bam_path = outdir / "possorted_genome_bam.bam"
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for query_name, cb, ub, cr, ur, mapq, unmapped, has_tx in reads:
            segment = pysam.AlignedSegment(out.header)
            segment.query_name = query_name
            segment.query_sequence = "ACGTACGTAC"
            segment.query_qualities = pysam.qualitystring_to_array("IIIIIIIIII")
            tags = [("CB", cb), ("UB", ub), ("CR", cr), ("UR", ur)]
            if has_tx:
                tags.append(("TX", "ENST00000000001,+100,10M"))
            if unmapped:
                segment.flag = 4
                segment.reference_id = -1
                segment.reference_start = -1
                segment.mapping_quality = 0
                segment.cigarstring = None
            else:
                segment.flag = 0
                segment.reference_id = 0
                segment.reference_start = 100
                segment.mapping_quality = mapq
                segment.cigarstring = "10M"
            segment.set_tags(tags)
            out.write(segment)

    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))

    expected = {
        "reads": [5, 1, 2, 0],
        "conf_mapped": [2, 1, 1, 0],
        "nonconf_mapped_reads": [1, 0, 1, 0],
        "unmapped_reads": [1, 0, 0, 0],
        "barcode_corrected_reads": [1, 0, 0, 0],
        "umi_corrected_reads": [1, 0, 0, 0],
        "barcode": [enc(barcode_0), enc(barcode_0), enc(barcode_1), enc(barcode_1)],
    }
    (outdir / "expected.json").write_text(json.dumps(expected, indent=2))
    print("Generated test data in", outdir)
    print("Expected:", json.dumps(expected))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
