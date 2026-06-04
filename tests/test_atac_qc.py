import pandas as pd
import pytest

from cellranger_qc.atac_qc import QCMetricsEngine, add_metadata


def test_add_metadata_accepts_bc_column(tmp_path):
    qc_metrics = pd.DataFrame(
        {
            "barcode": ["AAAC-1", "TTGG-1"],
            "atac_tss_enrichment": [6.0, 8.0],
        }
    )
    pd.DataFrame(
        {
            "bc": ["AAAC", "TTGG"],
            "cell_member": ["load_AAAC", "load_TTGG"],
        }
    ).to_csv(tmp_path / "samp.dat_AR123.csv", index=False)

    result = add_metadata(qc_metrics, str(tmp_path), "AR123")

    assert result["cell_member"].tolist() == ["load_AAAC", "load_TTGG"]
    assert "raw_barcode" not in result.columns


def test_add_metadata_accepts_barcodes_column(tmp_path):
    qc_metrics = pd.DataFrame({"barcode": ["AAAC-1"]})
    pd.DataFrame(
        {
            "barcodes": ["AAAC"],
            "cell_member": ["load_AAAC"],
        }
    ).to_csv(tmp_path / "samp.dat_AR123.csv", index=False)

    result = add_metadata(qc_metrics, str(tmp_path), "AR123")

    assert result.loc[0, "cell_member"] == "load_AAAC"


def test_add_metadata_requires_barcode_identifier(tmp_path):
    qc_metrics = pd.DataFrame({"barcode": ["AAAC-1"]})
    pd.DataFrame({"cell_member": ["load_AAAC"]}).to_csv(
        tmp_path / "samp.dat_AR123.csv", index=False
    )

    with pytest.raises(ValueError, match="must include either"):
        add_metadata(qc_metrics, str(tmp_path), "AR123")


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

    result = engine.compute_qc_metrics("load", nucleosome_df, insertion_counts_df)

    assert result.loc[0, "atac_nucleosome_ratio"] == 2
    assert result.loc[0, "atac_n_nucleosome_free_frags"] == 2
