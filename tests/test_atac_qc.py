import pandas as pd
import pytest

from cellranger_qc.atac_qc import add_metadata


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
