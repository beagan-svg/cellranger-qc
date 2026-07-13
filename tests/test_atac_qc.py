from contextlib import nullcontext as does_not_raise
from pathlib import Path

import pandas as pd
import pytest

from cellranger_qc.atac_qc import QCMetricsEngine, save_to_csv, validate_tabix_index


def test__save_to_csv__writes_standalone_output_name(tmp_path: Path):
    qc_metrics_df = pd.DataFrame({"barcode": ["AAAC-1"], "atac_tss_enrichment": [6.0]})
    expected = pd.DataFrame({"barcode": ["AAAC-1"], "atac_tss_enrichment": [6.0]})

    save_to_csv(qc_metrics_df, str(tmp_path))

    output_path = tmp_path / "atac_qc.csv"
    assert output_path.exists()
    assert expected.equals(pd.read_csv(output_path))


@pytest.mark.parametrize(
    "index_exists, raise_expectation",
    [
        pytest.param(
            # index_exists
            True,
            # raise_expectation
            does_not_raise(),
            id="Sibling tabix index is present",
        ),
        pytest.param(
            # index_exists
            False,
            # raise_expectation
            pytest.raises(FileNotFoundError, match="tabix-indexed"),
            id="Sibling tabix index is missing",
        ),
    ],
)
def test__validate_tabix_index__validates_sibling_index(
    tmp_path: Path,
    index_exists: bool,
    raise_expectation,
):
    fragments_path = tmp_path / "fragments.tsv.gz"
    fragments_path.touch()
    if index_exists:
        (tmp_path / "fragments.tsv.gz.tbi").touch()

    with raise_expectation:
        validate_tabix_index(str(fragments_path))


@pytest.mark.parametrize(
    (
        "nucleosome_df, insertion_counts_df, min_tss, min_frags_per_cell, "
        "raise_expectation, expected"
    ),
    [
        pytest.param(
            # nucleosome_df
            pd.DataFrame(
                {
                    "barcode": ["AAAC-1"],
                    "n_frags": [10],
                    "n_nucleosome_free_frags": [2],
                    "n_mono_frags": [3],
                    "n_di_frags": [1],
                    "n_multi_frags": [0],
                }
            ),
            # insertion_counts_df
            pd.DataFrame(
                {
                    "barcode": ["AAAC-1"],
                    "window": [101],
                    "flank": [100],
                    "promoter": [4],
                }
            ),
            # min_tss
            1,
            # min_frags_per_cell
            1,
            # raise_expectation
            does_not_raise(),
            # expected
            {"atac_nucleosome_ratio": 2.0, "atac_n_nucleosome_free_frags": 2},
            id="Nucleosome ratio uses the nucleosome-free fragment count",
        ),
        pytest.param(
            # nucleosome_df
            pd.DataFrame(
                {
                    "barcode": ["AAAC-1"],
                    "n_frags": [10],
                    "n_nucleosome_free_frags": [2],
                    "n_mono_frags": [1],
                    "n_di_frags": [0],
                    "n_multi_frags": [0],
                }
            ),
            # insertion_counts_df
            pd.DataFrame(
                {
                    "barcode": ["AAAC-1"],
                    "window": [10],
                    "flank": [100],
                    "promoter": [2],
                }
            ),
            # min_tss
            5,
            # min_frags_per_cell
            100,
            # raise_expectation
            pytest.raises(ValueError, match="No cells passed"),
            # expected
            None,
            id="No cells meet the configured filters",
        ),
    ],
)
def test__QCMetricsEngine__compute_qc_metrics(
    nucleosome_df: pd.DataFrame,
    insertion_counts_df: pd.DataFrame,
    min_tss: int,
    min_frags_per_cell: int,
    raise_expectation,
    expected: dict[str, float | int] | None,
):
    engine = object.__new__(QCMetricsEngine)
    engine.window = 101
    engine.norm = 100
    engine.min_norm = 0.2
    engine.min_tss = min_tss
    engine.min_frags_per_cell = min_frags_per_cell
    engine.max_frags_per_cell = 1000

    with raise_expectation:
        obt = engine.compute_qc_metrics(nucleosome_df, insertion_counts_df)

    if expected is not None:
        assert expected == {column: obt.loc[0, column] for column in expected}


def test__QCMetricsEngine__create_tss_regions__deduplicates_transcripts():
    engine = object.__new__(QCMetricsEngine)
    engine.window = 101
    engine.flank = 2000
    engine.norm = 100
    transcript_df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [5000, 5000],
            "end": [9000, 9000],
            "strand": ["+", "+"],
        }
    )
    expected = pd.DataFrame(
        [
            {"chrom": "chr1", "start": 3000, "end": 3099, "type": "flank"},
            {"chrom": "chr1", "start": 4950, "end": 5050, "type": "window"},
            {"chrom": "chr1", "start": 6901, "end": 7000, "type": "flank"},
        ]
    )

    obt = engine.create_tss_regions(transcript_df)

    assert expected.equals(obt)


def test__QCMetricsEngine__create_promoter_regions__uses_gene_strand():
    engine = object.__new__(QCMetricsEngine)
    genes_df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [5000, 5000],
            "end": [9000, 9000],
            "strand": ["+", "-"],
        }
    )
    expected = pd.DataFrame(
        [
            {"chrom": "chr1", "start": 3000, "end": 5100, "type": "promoter"},
            {"chrom": "chr1", "start": 8900, "end": 11000, "type": "promoter"},
        ]
    )

    obt = engine.create_promoter_regions(genes_df, region_span=(2000, 100))

    assert expected.equals(obt)


@pytest.mark.parametrize(
    "skip_chr_m, expected",
    [
        pytest.param(
            # skip_chr_m
            True,
            # expected
            ["chr1"],
            id="Mitochondrial chromosome is excluded",
        ),
        pytest.param(
            # skip_chr_m
            False,
            # expected
            ["chr1", "chrM"],
            id="Mitochondrial chromosome is retained",
        ),
    ],
)
def test__QCMetricsEngine__exclude_scaffold_chromosomes(
    skip_chr_m: bool,
    expected: list[str],
):
    engine = object.__new__(QCMetricsEngine)
    engine.skip_chr_m = skip_chr_m
    annotation_df = pd.DataFrame({"chrom": ["chr1", "chrM", "chrUn_1", "NW_123"]})

    obt = engine._exclude_scaffold_chromosomes(annotation_df)

    assert expected == obt["chrom"].tolist()
