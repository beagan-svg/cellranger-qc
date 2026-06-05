import gzip
import json
import shutil

import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse

from cellranger_qc.post_align_qc import (
    LoadedLibrary,
    get_cell_samp_dat,
    get_total_reads,
    load_data,
    write_summary_stats,
)


def test_get_cell_samp_dat_builds_expected_columns(monkeypatch, tmp_path):
    row_indices = np.concatenate([np.arange(1600), np.arange(1200)])
    column_indices = np.concatenate([np.zeros(1600, dtype=int), np.ones(1200, dtype=int)])
    count_matrix = sparse.csc_array(
        (np.ones(2800), (row_indices, column_indices)),
        shape=(1600, 2),
    )
    library = LoadedLibrary(
        count_matrix=count_matrix,
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
        umi_counts=np.asarray(count_matrix.sum(axis=0)).ravel(),
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


def test_get_total_reads_uses_per_barcode_metrics_when_available(tmp_path):
    outs_dir = tmp_path / "outs"
    outs_dir.mkdir()
    pd.DataFrame(
        {
            "gex_barcode": ["AAAC-1", "TTGG-1", "CCAA-2"],
            "is_cell": [1, 0, 1],
            "gex_raw_reads": [100, 200, 300],
        }
    ).to_csv(outs_dir / "per_barcode_metrics.csv", index=False)

    total_reads_df = get_total_reads(outs_dir)

    assert total_reads_df.to_dict("records") == [
        {"bc": "AAAC", "total_reads": 100},
        {"bc": "CCAA-2", "total_reads": 300},
    ]


def test_write_summary_stats_extracts_library_metrics(tmp_path):
    run_dir = tmp_path / "run"
    outs_dir = run_dir / "outs"
    outs_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "alignment_method": ["CELL_RANGER_COUNT"],
            "library_prep_method": ["GEX"],
            "Mean Reads per Cell": [1000],
        }
    ).to_csv(outs_dir / "summary.csv", index=False)

    web_summary = {"summary": {"diagnostics": {"tso_frac": 0.12}}}
    web_summary_line_list = [""] * 13
    web_summary_line_list[12] = "x" * 12 + json.dumps(web_summary)
    (outs_dir / "web_summary.html").write_text("\n".join(web_summary_line_list))

    samp_dat = pd.DataFrame(
        {
            "exclude": ["No", "No", "YES"],
            "exclude2": ["No", "YES", "YES"],
            "total_reads": [1000, 2000, 500],
            "gene_counts_0": [1500, 1000, 100],
        }
    )
    library_row = {
        "library_prep": "lib",
        "ar_dir": str(run_dir),
        "alignment_method": "CELL_RANGER_COUNT",
        "library_prep_method": "GEX",
        "expc_cell_capture": 4,
    }

    summary_df = write_summary_stats(samp_dat, library_row)

    assert summary_df.loc[0, "keeper_cells"] == 1
    assert summary_df.loc[0, "percent_keeper"] == 1 / 3
    assert summary_df.loc[0, "percent_doublet"] == 1 / 3
    assert summary_df.loc[0, "percent_usable"] == 0.25
    assert summary_df.loc[0, "tso_frac"] == 0.12
    assert summary_df.loc[0, "mean_reads_per_cell"] == 1000


def test_load_data_filters_multiome_gene_expression_and_disambiguates_duplicates(tmp_path):
    matrix_dir = tmp_path / "run" / "outs" / "filtered_feature_bc_matrix"
    matrix_dir.mkdir(parents=True)

    matrix_path = matrix_dir / "matrix.mtx"
    scipy.io.mmwrite(matrix_path, sparse.coo_array([[1, 2], [3, 4], [5, 6]]))
    with matrix_path.open("rb") as source, gzip.open(matrix_dir / "matrix.mtx.gz", "wb") as target:
        shutil.copyfileobj(source, target)
    matrix_path.unlink()

    pd.DataFrame(
        [
            ["gene-a", "GeneA", "Gene Expression"],
            ["peak-a", "PeakA", "Peaks"],
            ["gene-b", "GeneA", "Gene Expression"],
        ]
    ).to_csv(matrix_dir / "features.tsv.gz", sep="\t", header=False, index=False)
    pd.Series(["AAAC-1", "TTGG-1"]).to_csv(matrix_dir / "barcodes.tsv.gz", header=False, index=False)

    loaded_library = load_data(
        {
            "ar_dir": str(tmp_path / "run"),
            "alignment_method": "CELL_RANGER_MULTI",
            "library_prep": "lib",
            "ar_id": "AR123",
        }
    )

    assert loaded_library.count_matrix.shape == (2, 2)
    assert loaded_library.barcode_list.tolist() == ["AAAC", "TTGG"]
    assert loaded_library.gene_names.tolist() == ["GeneA", "GeneA gene-b"]
    assert loaded_library.sample_id.tolist() == ["AAAC-lib-AR123", "TTGG-lib-AR123"]
