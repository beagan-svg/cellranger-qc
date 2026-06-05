import pandas as pd
import pytest

from cellranger_qc.atac_qc import QCMetricsEngine, save_to_csv, validate_tabix_index


def test_save_to_csv_writes_standalone_output_name(tmp_path):
    qc_metrics = pd.DataFrame({"barcode": ["AAAC-1"], "atac_tss_enrichment": [6.0]})

    save_to_csv(qc_metrics, str(tmp_path))

    output_path = tmp_path / "atac_qc.csv"
    assert output_path.exists()
    assert pd.read_csv(output_path).to_dict("records") == [
        {"barcode": "AAAC-1", "atac_tss_enrichment": 6.0}
    ]


def test_validate_tabix_index_accepts_sibling_tbi(tmp_path):
    fragments_path = tmp_path / "fragments.tsv.gz"
    fragments_path.touch()
    (tmp_path / "fragments.tsv.gz.tbi").touch()

    validate_tabix_index(str(fragments_path))


def test_validate_tabix_index_requires_sibling_tbi(tmp_path):
    fragments_path = tmp_path / "fragments.tsv.gz"
    fragments_path.touch()

    with pytest.raises(FileNotFoundError, match="tabix-indexed"):
        validate_tabix_index(str(fragments_path))


def test_compute_qc_metrics_uses_nucleosome_free_denominator():
    engine = object.__new__(QCMetricsEngine)
    engine.window = 101
    engine.norm = 100
    engine.min_norm = 0.2
    engine.min_tss = 1
    engine.min_frags_per_cell = 1
    engine.max_frags_per_cell = 1000

    nucleosome_df = pd.DataFrame(
        {
            "barcode": ["AAAC-1"],
            "n_frags": [10],
            "n_nucleosome_free_frags": [2],
            "n_mono_frags": [3],
            "n_di_frags": [1],
            "n_multi_frags": [0],
        }
    )
    insertion_counts_df = pd.DataFrame(
        {
            "barcode": ["AAAC-1"],
            "window": [101],
            "flank": [100],
            "promoter": [4],
        }
    )

    result = engine.compute_qc_metrics(nucleosome_df, insertion_counts_df)

    assert result.loc[0, "atac_nucleosome_ratio"] == 2
    assert result.loc[0, "atac_n_nucleosome_free_frags"] == 2
