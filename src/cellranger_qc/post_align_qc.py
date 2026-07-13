import argparse
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse, stats
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from cellranger_qc import __version__

logger = logging.getLogger(__name__)

MULTIOME_ALIGNMENT_METHODS = {"CELL_RANGER_ARCV1", "ARC-RSEQ", "CELL_RANGER_MULTI"}
BARCODE_SUFFIX_PATTERN = r"-1$"
GENE_COUNT_THRESHOLDS = (0, 1, 4, 8, 16, 32, 64)
MAX_VARIABLE_GENES = 5000


def calculate_doublets(
    count_matrix: sparse.csc_array,
    gene_mask_indices: np.ndarray,
    sample_id: np.ndarray,
    doublet_out_path: Path | str | None = None,
    proportion_artificial: float = 0.20,
) -> pd.Series:
    """Score cells by their proximity to synthetic doublets in PCA space."""
    logger.info("Creating synthetic doublets")
    rng = np.random.RandomState(1)
    real_cell_count = count_matrix.shape[1]
    k = max(1, round(min(100, real_cell_count * 0.01)))

    doublet_count = round(real_cell_count / (1 - proportion_artificial) - real_cell_count)
    first_cell_indices = rng.choice(real_cell_count, doublet_count, replace=True)
    second_cell_indices = rng.choice(real_cell_count, doublet_count, replace=True)
    doublets = count_matrix[:, first_cell_indices] + count_matrix[:, second_cell_indices]
    data_with_doublets = sparse.hstack([count_matrix, doublets], format="csc")
    logger.info(f"Created {doublet_count} synthetic doublets")

    # log-CPM normalize per cell to account for sequencing depth
    cell_sums = np.asarray(data_with_doublets.sum(axis=0)).ravel()
    size_factors = np.divide(
        1e6,
        cell_sums,
        out=np.zeros(cell_sums.shape, dtype=float),
        where=cell_sums != 0,
    )
    norm_data = (data_with_doublets @ sparse.diags(size_factors, format="csc")).tocsc()
    norm_data.data = np.log2(norm_data.data + 1)

    logger.info("Embedding real + artificial cells together in PCA space")
    pca_input = norm_data[gene_mask_indices, :].T
    fit_indices = np.arange(pca_input.shape[0])
    if real_cell_count > 10000:
        fit_indices = rng.choice(real_cell_count, 10000, replace=False)
    pca = PCA(svd_solver="full", random_state=1)
    pca.fit(pca_input[fit_indices].toarray())
    # Keep components capturing meaningful variance (positive z-score), capped at 50
    component_mask = stats.zscore(pca.explained_variance_) > 0
    component_mask[50:] = False
    reduced_data = pca.transform(pca_input.toarray())[:, component_mask]

    logger.info("Scoring cells by synthetic-neighbor enrichment")
    knn = NearestNeighbors(n_neighbors=k).fit(reduced_data)
    neighbor_distances, neighbor_indices = knn.kneighbors(reduced_data)
    artificial_knn = NearestNeighbors(n_neighbors=10).fit(reduced_data[:real_cell_count])
    artificial_distances, _ = artificial_knn.kneighbors(reduced_data[real_cell_count:])
    distance_threshold = artificial_distances.mean() + 1.64 * artificial_distances.std(ddof=1)

    doublet_neighborhood = (neighbor_indices >= real_cell_count) & (
        neighbor_distances < distance_threshold
    )
    half_k = int(np.ceil(k / 2))
    doublet_score = np.maximum(
        doublet_neighborhood.mean(axis=1),
        doublet_neighborhood[:, :half_k].mean(axis=1),
    )[:real_cell_count]

    if doublet_out_path is not None:
        plt.figure()
        plt.hist(doublet_score, bins=100, density=True)
        plt.title("Doublet score density")
        plt.xlabel("Doublet score")
        plt.ylabel("Density")
        plt.tight_layout()
        plt.savefig(doublet_out_path)
        plt.close()

    return pd.Series(doublet_score, index=sample_id)


def get_total_reads(outs_dir: Path) -> pd.DataFrame:
    """Read total per-cell counts from Cell Ranger barcode metrics or molecule info."""
    per_barcode_metrics_path = outs_dir / "per_barcode_metrics.csv"
    if per_barcode_metrics_path.exists():
        barcode_metrics = pd.read_csv(per_barcode_metrics_path)
        barcode_metrics = barcode_metrics.loc[barcode_metrics["is_cell"] == 1]
        return pd.DataFrame(
            {
                "bc": barcode_metrics["gex_barcode"].str.replace(
                    BARCODE_SUFFIX_PATTERN, "", regex=True
                ),
                "total_reads": barcode_metrics["gex_raw_reads"],
            }
        )

    molecule_info_path = next(outs_dir.glob("*molecule_info_new.h5"))
    with h5py.File(molecule_info_path, "r") as h5_file:
        barcodes = h5_file["barcodes"][...].astype(str)
        barcode_indices = h5_file["barcode_idx"][...]
        total_reads = (
            h5_file["reads"][...]
            + h5_file["unmapped_reads"][...]
            + h5_file["nonconf_mapped_reads"][...]
        )
        pass_filter_indices = h5_file["barcode_info/pass_filter"][:, 0]

    reads = pd.DataFrame({"barcode_indices": barcode_indices, "total_reads": total_reads})
    reads = reads.loc[reads["barcode_indices"].isin(pass_filter_indices)]
    reads = reads.groupby("barcode_indices", as_index=False)["total_reads"].sum()
    reads["bc"] = pd.Series(barcodes[reads["barcode_indices"].to_numpy()]).str.replace(
        BARCODE_SUFFIX_PATTERN, "", regex=True
    )
    return reads[["bc", "total_reads"]]


def get_cell_samp_dat(
    loaded_library: "LoadedLibrary",
    umi_counts: np.ndarray,
    library_row: dict[str, Any],
    out_dir: Path | str,
) -> pd.DataFrame:
    """Build the per-cell QC table for one library."""
    samp_dat = pd.DataFrame(
        {
            "sample_id": loaded_library.sample_id,
            "bc": loaded_library.barcode_list,
            "umi_counts": umi_counts,
            "library_prep": loaded_library.library_prep,
        }
    )
    for gene_threshold in GENE_COUNT_THRESHOLDS:
        samp_dat[f"gene_counts_{gene_threshold}"] = (
            loaded_library.count_matrix > gene_threshold
        ).sum(axis=0)

    count_matrix = loaded_library.count_matrix.copy()
    count_matrix.data = np.log2(count_matrix.data + 1)

    # Variance per gene across all cells, used to pick the top 5,000 variable genes
    gene_means = count_matrix.mean(axis=1)
    gene_variance = count_matrix.power(2).mean(axis=1) - gene_means**2
    top_gene_indices = np.argsort(gene_variance)[::-1][:MAX_VARIABLE_GENES]

    doublet_out_path = Path(out_dir) / f"{loaded_library.library_prep}.doubscore.pdf"
    doublet_score = calculate_doublets(
        count_matrix, top_gene_indices, loaded_library.sample_id, doublet_out_path
    )
    samp_dat["doublet_score"] = samp_dat["sample_id"].map(doublet_score)

    gene_count_threshold = 1500 if library_row["cell_prep_type"] == "Cells" else 1000
    samp_dat["exclude"] = np.where(samp_dat["gene_counts_0"] < gene_count_threshold, "YES", "No")
    samp_dat["exclude2"] = np.where(
        (samp_dat["exclude"] == "YES") | (samp_dat["doublet_score"] > 0.3), "YES", "No"
    )

    outs_dir = Path(library_row["cellranger_run_dir"]) / "outs"
    samp_dat = samp_dat.merge(get_total_reads(outs_dir), on="bc")
    samp_dat["cell_member"] = str(library_row["load_name"]) + "_" + samp_dat["bc"]
    return samp_dat


def write_summary_stats(samp_dat: pd.DataFrame, library_row: dict[str, Any]) -> pd.DataFrame:
    """Combine Cell Ranger metrics with library-level keeper statistics."""
    keepers = samp_dat.loc[samp_dat["exclude"] == "No"]
    keeper_cells = int((keepers["exclude2"] == "No").sum())

    outs_dir = Path(library_row["cellranger_run_dir"]) / "outs"
    alignment_metrics = pd.read_csv(next(outs_dir.glob("*summary.csv")))
    alignment_metrics.columns = alignment_metrics.columns.str.replace(
        r"[^0-9A-Za-z_]+", "_", regex=True
    ).str.lower()
    library_summary = pd.DataFrame([library_row])
    library_summary = library_summary[
        ["library_prep"]
        + [column for column in library_summary.columns if column != "library_prep"]
    ]

    ocs_summary = pd.concat(
        [
            library_summary.reset_index(drop=True),
            alignment_metrics.reset_index(drop=True),
            pd.DataFrame(
                [
                    {
                        "keeper_mean": keepers["total_reads"].mean(),
                        "keeper_median_genes": keepers["gene_counts_0"].median(),
                        "keeper_cells": keeper_cells,
                        "percent_keeper": keeper_cells / len(samp_dat),
                        "percent_doublet": (len(keepers) - keeper_cells) / len(samp_dat),
                        "percent_usable": keeper_cells / library_row["expc_cell_capture"],
                    }
                ]
            ),
        ],
        axis=1,
    )
    ocs_summary = ocs_summary.drop(columns=["alignment_method", "library_prep_method"])

    web_summary_lines = (outs_dir / "web_summary.html").read_text().splitlines()
    if library_row["alignment_method"] in {"ARC-RSEQ", "CELL_RANGER_MULTI"}:
        web_summary = json.loads(web_summary_lines[222].strip()[12:])
        ocs_summary["tso_frac"] = (
            float(web_summary["gex_sequencing_table"]["rows"][8][1].replace("%", "")) / 100
        )
    else:
        web_summary = json.loads(web_summary_lines[12].strip()[12:])
        ocs_summary["tso_frac"] = web_summary["summary"]["diagnostics"]["tso_frac"]

    ocs_summary["pass_fail"] = "pass"
    ocs_summary.columns = ocs_summary.columns.str.replace(
        r"[^0-9A-Za-z_]+", "_", regex=True
    ).str.lower()
    return ocs_summary


def extract_intron_exon_matrices(molecule_info_path: Path) -> dict[str, sparse.csc_array]:
    """Build pass-filtered intron and exon matrices from `molecule_info.h5`."""
    with h5py.File(molecule_info_path, "r") as h5_file:
        barcodes = np.char.replace(h5_file["barcodes"][...].astype(str), "-1", "")
        molecule_barcode_indices = h5_file["barcode_idx"][...]
        molecule_feature_indices = h5_file["feature_idx"][...]
        molecule_umi_type = h5_file["umi_type"][...]

        feature_names = h5_file["features/name"][...].astype(str).tolist()
        n_features = len(feature_names)

        pass_filter_table = h5_file["barcode_info/pass_filter"][...]
        pass_filter_barcode_indices = pass_filter_table[:, 0].astype(int)

    filtered_barcodes = barcodes[pass_filter_barcode_indices]

    barcode_index_to_column = np.full(len(barcodes), -1, dtype=np.int32)
    barcode_index_to_column[pass_filter_barcode_indices] = np.arange(
        len(pass_filter_barcode_indices), dtype=np.int32
    )

    is_pass_filter_molecule = barcode_index_to_column[molecule_barcode_indices] >= 0

    matrix_shape = (n_features, len(filtered_barcodes))

    def build_umi_matrix(umi_value: int) -> sparse.csc_array:
        """Build the matrix for one Cell Ranger UMI type."""
        is_matching_molecule = is_pass_filter_molecule & (molecule_umi_type == umi_value)
        return sparse.coo_array(
            (
                np.ones(is_matching_molecule.sum()),
                (
                    molecule_feature_indices[is_matching_molecule],
                    barcode_index_to_column[molecule_barcode_indices[is_matching_molecule]],
                ),
            ),
            shape=matrix_shape,
        ).tocsc()

    return {
        "exons": build_umi_matrix(1),
        "introns": build_umi_matrix(0),
    }


def generate_intron_exon(
    cellranger_run_dir: Path | str, output_prefix: str, out_dir: Path | str
) -> None:
    """Write intron and exon matrices for one Cell Ranger run."""
    outs_dir = Path(cellranger_run_dir) / "outs"
    molecule_info_path = next(outs_dir.glob("*molecule_info.h5"), None)
    if molecule_info_path is None:
        raise FileNotFoundError(
            f"No *molecule_info.h5 in {outs_dir}\n"
            + "\n".join(f"  {path.name}" for path in outs_dir.iterdir())
        )

    intron_exon_matrices = extract_intron_exon_matrices(molecule_info_path)
    logger.info("Saving exon/intron matrices")
    scipy.io.mmwrite(
        Path(out_dir) / "matrix" / f"intron_{output_prefix}.mtx", intron_exon_matrices["introns"]
    )
    scipy.io.mmwrite(
        Path(out_dir) / "matrix" / f"exon_{output_prefix}.mtx", intron_exon_matrices["exons"]
    )


@dataclass
class LoadedLibrary:
    """Filtered count data and identifiers for one Cell Ranger library."""

    count_matrix: sparse.csc_array
    gene_df: pd.DataFrame
    barcode_list: np.ndarray
    sample_id: np.ndarray
    gene_names: np.ndarray
    library_prep: str


def load_data(library_row: dict[str, Any]) -> LoadedLibrary:
    """Load and normalize filtered feature-barcode data for one library."""
    matrix_dir = Path(library_row["cellranger_run_dir"]) / "outs" / "filtered_feature_bc_matrix"
    count_matrix = sparse.csc_array(scipy.io.mmread(matrix_dir / "matrix.mtx.gz"))
    gene_df = pd.read_csv(matrix_dir / "features.tsv.gz", sep="\t", header=None)
    barcode_list = pd.read_csv(matrix_dir / "barcodes.tsv.gz", header=None)[0].str.replace(
        BARCODE_SUFFIX_PATTERN, "", regex=True
    )

    if library_row["alignment_method"] in MULTIOME_ALIGNMENT_METHODS:
        is_gene_expression = gene_df[2].eq("Gene Expression").to_numpy()
        gene_df = gene_df.loc[is_gene_expression].reset_index(drop=True)
        count_matrix = count_matrix[is_gene_expression, :]

    gene_names_series = gene_df[1].astype(str)
    is_duplicate_gene_name = gene_names_series.duplicated()
    gene_names_series.loc[is_duplicate_gene_name] = (
        gene_df.loc[is_duplicate_gene_name, 1].astype(str)
        + " "
        + gene_df.loc[is_duplicate_gene_name, 0].astype(str)
    )

    sample_id = barcode_list + "-" + str(library_row["library_prep"])

    return LoadedLibrary(
        count_matrix=count_matrix.tocsc(),
        gene_df=gene_df,
        barcode_list=barcode_list.to_numpy(),
        sample_id=sample_id.to_numpy(),
        gene_names=gene_names_series.to_numpy(),
        library_prep=library_row["library_prep"],
    )


def run_rseq_qc(libs: pd.DataFrame, out_dir: Path | str, num_cores: int = 16) -> None:
    """Run GEX or multiome QC for every library in the manifest."""
    out_dir = Path(out_dir)
    (out_dir / "matrix").mkdir(parents=True, exist_ok=True)

    logger.info(f"out_dir: {out_dir}/")
    logger.info(f"number of libraries =>  {len(libs)}")
    logger.info(f"number of cores =>  {num_cores}")

    for _, library in libs.iterrows():
        library_row = library.to_dict()

        loaded_library = load_data(library_row)
        logger.info(
            f"{loaded_library.count_matrix.shape[0]} genes x "
            f"{loaded_library.count_matrix.shape[1]} cells"
        )
        output_prefix = str(library_row["library_prep"])
        generate_intron_exon(library_row["cellranger_run_dir"], output_prefix, out_dir)
        logger.info("Saving count matrix")
        scipy.io.mmwrite(
            out_dir / "matrix" / f"count_{output_prefix}.mtx",
            loaded_library.count_matrix,
        )

        umi_counts = loaded_library.count_matrix.sum(axis=0)
        samp_dat = get_cell_samp_dat(loaded_library, umi_counts, library_row, out_dir)
        samp_dat.to_csv(out_dir / f"samp_dat_{output_prefix}.csv", index=False)

        ocs_summary = write_summary_stats(samp_dat, library_row)
        ocs_summary.to_csv(out_dir / "ocs_summary.csv", index=False)


def main() -> None:
    """Run the GEX QC command-line interface."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    parser = argparse.ArgumentParser()
    parser.add_argument("--libs", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--num-cores", type=int, default=16)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()

    libs = pd.read_csv(args.libs)
    run_rseq_qc(libs, args.out_dir, args.num_cores)


if __name__ == "__main__":
    main()
