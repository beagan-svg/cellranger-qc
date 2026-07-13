import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.io
from scipy import sparse

from cellranger_qc.post_align_qc import (
    LoadedLibrary,
    get_cell_samp_dat,
    get_total_reads,
    load_data,
    write_summary_stats,
)


def test__get_cell_samp_dat__builds_expected_exclusion_columns(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    row_indices = np.concatenate([np.arange(1600), np.arange(1200)])
    column_indices = np.concatenate([np.zeros(1600, dtype=int), np.ones(1200, dtype=int)])
    count_matrix = sparse.csc_array(
        (np.ones(2800), (row_indices, column_indices)),
        shape=(1600, 2),
    )
    loaded_library = LoadedLibrary(
        count_matrix=count_matrix,
        gene_df=pd.DataFrame(),
        barcode_list=np.array(["AAAC", "TTGG"]),
        sample_id=np.array(["AAAC-lib", "TTGG-lib"]),
        gene_names=np.array(["GeneA", "GeneB"]),
        library_prep="lib",
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
    expected = pd.DataFrame(
        {
            "exclude": ["No", "YES"],
            "exclude2": ["No", "YES"],
            "cell_member": ["load_AAAC", "load_TTGG"],
        }
    )

    obt = get_cell_samp_dat(
        loaded_library,
        umi_counts=np.asarray(count_matrix.sum(axis=0)).ravel(),
        library_row={
            "cell_prep_type": "Cells",
            "cellranger_run_dir": str(tmp_path),
            "load_name": "load",
        },
        out_dir=tmp_path,
    )

    assert expected.equals(obt[list(expected.columns)])


def test__get_total_reads__uses_per_barcode_metrics_when_available(tmp_path: Path):
    outs_dir = tmp_path / "outs"
    outs_dir.mkdir()
    pd.DataFrame(
        {
            "gex_barcode": ["AAAC-1", "TTGG-1", "CCAA-2"],
            "is_cell": [1, 0, 1],
            "gex_raw_reads": [100, 200, 300],
        }
    ).to_csv(outs_dir / "per_barcode_metrics.csv", index=False)
    expected = pd.DataFrame(
        [
            {"bc": "AAAC", "total_reads": 100},
            {"bc": "CCAA-2", "total_reads": 300},
        ]
    )

    obt = get_total_reads(outs_dir).reset_index(drop=True)

    assert expected.equals(obt)


def test__write_summary_stats__extracts_library_metrics(tmp_path: Path):
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
    web_summary_lines = [""] * 13
    web_summary_lines[12] = "x" * 12 + json.dumps(web_summary)
    (outs_dir / "web_summary.html").write_text("\n".join(web_summary_lines))

    samp_dat_df = pd.DataFrame(
        {
            "exclude": ["No", "No", "YES"],
            "exclude2": ["No", "YES", "YES"],
            "total_reads": [1000, 2000, 500],
            "gene_counts_0": [1500, 1000, 100],
        }
    )
    library_row = {
        "library_prep": "lib",
        "cellranger_run_dir": str(run_dir),
        "alignment_method": "CELL_RANGER_COUNT",
        "library_prep_method": "GEX",
        "expc_cell_capture": 4,
    }
    expected = {
        "keeper_cells": 1,
        "percent_keeper": 1 / 3,
        "percent_doublet": 1 / 3,
        "percent_usable": 0.25,
        "tso_frac": 0.12,
        "mean_reads_per_cell": 1000,
    }

    obt = write_summary_stats(samp_dat_df, library_row)

    assert expected == {column: obt.loc[0, column] for column in expected}


def test__load_data__filters_multiome_features_and_disambiguates_gene_names(tmp_path: Path):
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
    pd.Series(["AAAC-1", "TTGG-1"]).to_csv(
        matrix_dir / "barcodes.tsv.gz", header=False, index=False
    )

    obt = load_data(
        {
            "cellranger_run_dir": str(tmp_path / "run"),
            "alignment_method": "CELL_RANGER_MULTI",
            "library_prep": "lib",
        }
    )

    assert obt.count_matrix.shape == (2, 2)
    assert obt.barcode_list.tolist() == ["AAAC", "TTGG"]
    assert obt.gene_names.tolist() == ["GeneA", "GeneA gene-b"]
    assert obt.sample_id.tolist() == ["AAAC-lib", "TTGG-lib"]
