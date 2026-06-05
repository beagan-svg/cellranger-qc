import numpy as np
import pandas as pd

from cellranger_qc.post_align_qc import LoadedLibrary, get_cell_samp_dat


class FakeMatrix:
    def __init__(self, threshold_counts):
        self.threshold_counts = threshold_counts
        self.data = np.array([1.0])

    def __gt__(self, threshold):
        return FakeComparison(self.threshold_counts[threshold])

    def copy(self):
        copied_matrix = FakeMatrix(self.threshold_counts)
        copied_matrix.data = self.data.copy()
        return copied_matrix

    def mean(self, axis):
        return np.array([[0.0], [1.0], [2.0]])

    def power(self, value):
        return self


class FakeComparison:
    def __init__(self, values):
        self.values = np.asarray(values)

    def sum(self, axis):
        return self.values


def test_get_cell_samp_dat_builds_expected_columns(monkeypatch, tmp_path):
    matrix = FakeMatrix(
        {
            0: [1600, 1200],
            1: [1500, 1100],
            4: [1000, 900],
            8: [500, 400],
            16: [200, 100],
            32: [50, 25],
            64: [10, 5],
        }
    )
    library = LoadedLibrary(
        count_matrix=matrix,
        gene_df=pd.DataFrame(),
        barcode_list=np.array(["AAAC", "TTGG"]),
        sample_id=np.array(["AAAC-lib-AR123", "TTGG-lib-AR123"]),
        gene_names=np.array(["GeneA", "GeneB"]),
        library_prep="lib",
        ar_id="AR123",
    )

    monkeypatch.setattr(
        "cellranger_qc.post_align_qc.calculate_doublets",
        lambda count_matrix, gene_mask_indices, sample_id, doublet_out_path: pd.Series(
            [0.1, 0.4], index=sample_id
        ),
    )
    monkeypatch.setattr(
        "cellranger_qc.post_align_qc.get_total_reads",
        lambda outs_dir: pd.DataFrame({"bc": ["AAAC", "TTGG"], "total_reads": [10000, 8000]}),
    )

    result = get_cell_samp_dat(
        library,
        umi_counts=np.array([5000, 4000]),
        library_row={
            "cell_prep_type": "Cells",
            "ar_dir": str(tmp_path),
            "load_name": "load",
        },
        out_dir=tmp_path,
    )

    assert result["exclude"].tolist() == ["No", "YES"]
    assert result["exclude2"].tolist() == ["No", "YES"]
    assert result["cell_member"].tolist() == ["load_AAAC", "load_TTGG"]
