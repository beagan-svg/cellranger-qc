"""Calculate ATAC QC metrics from Cell Ranger fragments and annotations.

Fragment coordinates are BED-like. GTF coordinates are one-based and inclusive.
Each fragment contributes its two Tn5 insertion endpoints to overlapping regions.
"""

import argparse
import logging
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pysam
from ncls import NCLS

from cellranger_qc import __version__

MAX_INSERTION_WORKERS = 8

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logger = logging.getLogger(__name__)


class QCMetricsEngine:
    """Compute TSS enrichment, promoter ratio, and nucleosome metrics."""

    def __init__(
        self,
        annotation_file: str,
        atac_fragments_path: str,
        per_barcode_metrics_path: str,
        window: int = 101,
        flank: int = 2000,
        norm: int = 100,
        min_norm: float = 0.2,
        skip_chr_m: bool = True,
        min_tss: float = 1,
        min_frags_per_cell: int = 1000,
        max_frags_per_cell: int = 100000000,
        nuc_min_frags: int = 10,
        nuc_max_frags: int = 2000,
        nuc_len: int = 147,
    ):
        """Load cell barcodes and configure the ATAC QC thresholds."""
        self.annotation_file = annotation_file
        self.atac_fragments_path = atac_fragments_path
        per_barcode_metrics_df = pd.read_csv(
            per_barcode_metrics_path,
            usecols=["barcode", "is_cell"],
        )

        self.is_cell_bc = set(
            per_barcode_metrics_df[per_barcode_metrics_df["is_cell"] == 1]["barcode"]
        )
        total_barcodes = len(per_barcode_metrics_df)
        filtered_percent = len(self.is_cell_bc) / total_barcodes * 100 if total_barcodes else 0
        logger.info(
            "Filtered cells: %s out of %s total barcodes (%.1f%%)",
            len(self.is_cell_bc),
            total_barcodes,
            filtered_percent,
        )
        self.window = window
        self.flank = flank
        self.norm = norm
        self.min_norm = min_norm
        self.skip_chr_m = skip_chr_m
        self.min_tss = min_tss
        self.min_frags_per_cell = min_frags_per_cell
        self.max_frags_per_cell = max_frags_per_cell
        self.nuc_min_frags = nuc_min_frags
        self.nuc_max_frags = nuc_max_frags
        self.nuc_len = nuc_len
        self.annotation_regions = pd.DataFrame()
        self.valid_chromosomes = list()
        self.barcodes = sorted(self.is_cell_bc)
        self.barcode_to_row = {barcode: row for row, barcode in enumerate(self.barcodes)}
        logger.info(
            "ATAC QC parameters: window=%s, flank=%s, norm=%s, min_norm=%s, "
            "skip_chr_m=%s, min_tss=%s, min_frags_per_cell=%s, "
            "max_frags_per_cell=%s, nuc_min_frags=%s, nuc_max_frags=%s, nuc_len=%s",
            self.window,
            self.flank,
            self.norm,
            self.min_norm,
            self.skip_chr_m,
            self.min_tss,
            self.min_frags_per_cell,
            self.max_frags_per_cell,
            self.nuc_min_frags,
            self.nuc_max_frags,
            self.nuc_len,
        )

    def nucleosome_classification(self) -> pd.DataFrame:
        """Count nucleosome-free, mono-, di-, and multi-nucleosome fragments per cell."""
        logger.info("Fragments file path: %s", self.atac_fragments_path)

        file_size_bytes = Path(self.atac_fragments_path).stat().st_size
        file_size_mb = file_size_bytes / (1024 * 1024)
        file_size_gb = file_size_bytes / (1024 * 1024 * 1024)
        if file_size_gb >= 1:
            logger.info(
                "ATAC fragment file size: %.2f GB (%s bytes)",
                file_size_gb,
                f"{file_size_bytes:,}",
            )
        else:
            logger.info(
                "ATAC fragment file size: %.2f MB (%s bytes)",
                file_size_mb,
                f"{file_size_bytes:,}",
            )
        start_time = time.time()

        filtered_fragments_df = (
            pl.scan_csv(
                self.atac_fragments_path,
                has_header=False,
                separator="\t",
                comment_prefix="#",
            )
            .select(
                [
                    pl.col("column_1").cast(pl.Utf8).alias("chrom"),
                    pl.col("column_2").cast(pl.Int64).alias("start"),
                    pl.col("column_3").cast(pl.Int64).alias("end"),
                    pl.col("column_4").cast(pl.Utf8).alias("barcode"),
                ]
            )
            .drop_nulls(["chrom", "start", "end", "barcode"])
            .filter(pl.col("barcode").is_in(self.is_cell_bc))
            .filter(pl.col("chrom").is_in(self.valid_chromosomes))
            .with_columns(
                [
                    (pl.col("end") - pl.col("start")).alias("frag_size"),
                ]
            )
            .filter(
                (pl.col("frag_size") >= self.nuc_min_frags)
                & (pl.col("frag_size") <= self.nuc_max_frags)
            )
            .collect()
        )

        nucleosome_df = (
            filtered_fragments_df.group_by(["barcode"])
            .agg(
                [
                    pl.len().alias("n_frags"),
                    pl.col("frag_size")
                    .filter(pl.col("frag_size") < self.nuc_len)
                    .count()
                    .alias("n_nucleosome_free_frags"),
                    pl.col("frag_size")
                    .filter(
                        (pl.col("frag_size") >= self.nuc_len)
                        & (pl.col("frag_size") < 2 * self.nuc_len)
                    )
                    .count()
                    .alias("n_mono_frags"),
                    pl.col("frag_size")
                    .filter(
                        (pl.col("frag_size") >= 2 * self.nuc_len)
                        & (pl.col("frag_size") < 3 * self.nuc_len)
                    )
                    .count()
                    .alias("n_di_frags"),
                    pl.col("frag_size")
                    .filter(pl.col("frag_size") >= 3 * self.nuc_len)
                    .count()
                    .alias("n_multi_frags"),
                ]
            )
            .to_pandas()
        )

        logger.info("Total fragments: %s", f"{filtered_fragments_df.height:,}")
        nucleosome_time = time.time() - start_time
        logger.info(
            "Nucleosome counting processing time: %.2f sec (%.2f min)",
            nucleosome_time,
            nucleosome_time / 60,
        )

        return nucleosome_df

    def count_fragments_and_insertions(self) -> pd.DataFrame:
        """Count per-cell insertions in TSS windows, flanks, and promoters."""
        start_time = time.time()
        insertion_counts = np.zeros((len(self.barcodes), 3), dtype=np.int64)

        with pysam.TabixFile(self.atac_fragments_path) as tbx_file:
            available_contigs = set(tbx_file.contigs)

        chromosome_batches = list()
        for _ in range(MAX_INSERTION_WORKERS):
            chromosome_batches.append(list())
        batch_sizes = [0] * MAX_INSERTION_WORKERS
        for chromosome, chromosome_regions in self.annotation_regions.groupby("chrom", sort=False):
            if chromosome not in available_contigs:
                continue

            worker_index = batch_sizes.index(min(batch_sizes))
            chromosome_batches[worker_index].append((chromosome, chromosome_regions))
            batch_sizes[worker_index] += len(chromosome_regions)

        with ProcessPoolExecutor(max_workers=MAX_INSERTION_WORKERS) as executor:
            for chromosome_batch_counts in executor.map(
                self._count_insertions_for_chromosome_batch, chromosome_batches
            ):
                insertion_counts += chromosome_batch_counts

        insertion_counts_df = pd.DataFrame(
            {
                "barcode": self.barcodes,
                "window": insertion_counts[:, 0],
                "flank": insertion_counts[:, 1],
                "promoter": insertion_counts[:, 2],
            }
        )
        insertion_counts_df = insertion_counts_df[
            (insertion_counts_df["window"] > 0)
            | (insertion_counts_df["flank"] > 0)
            | (insertion_counts_df["promoter"] > 0)
        ]

        elapsed_time = time.time() - start_time
        logger.info("TSS and promoter insertion counting summary:")
        logger.info("  Window insertions: %s", f"{insertion_counts[:, 0].sum():,}")
        logger.info("  Flank insertions: %s", f"{insertion_counts[:, 1].sum():,}")
        logger.info(
            "  Promoter insertions (start + end points): %s",
            f"{insertion_counts[:, 2].sum():,}",
        )
        logger.info("  Processing time: %.1fs (%.1f min)", elapsed_time, elapsed_time / 60)
        return insertion_counts_df

    def _count_insertions_for_chromosome_batch(self, chromosome_batch) -> np.ndarray:
        """Count region-overlapping insertions for one worker's chromosome batch."""
        region_type_to_column = {"window": 0, "flank": 1, "promoter": 2}
        chromosome_batch_insertion_counts = np.zeros((len(self.barcodes), 3), dtype=np.int64)

        with pysam.TabixFile(self.atac_fragments_path) as tbx_file:
            for chromosome, chromosome_regions in chromosome_batch:
                chromosome_regions = chromosome_regions.sort_values(["start", "end"])
                region_starts = chromosome_regions["start"].to_numpy(dtype=np.int64)
                # NCLS uses half-open intervals.
                region_ends = chromosome_regions["end"].to_numpy(dtype=np.int64) + 1
                region_rows = np.arange(len(chromosome_regions), dtype=np.int64)
                region_columns = np.asarray(
                    [region_type_to_column[value] for value in chromosome_regions["type"]],
                    dtype=np.int8,
                )
                region_overlap_lookup = NCLS(region_starts, region_ends, region_rows)
                insertion_positions_list = list()
                barcode_indices_list = list()

                merged_start = region_starts[0]
                merged_end = region_ends[0]
                for region_start, region_end in zip(region_starts[1:], region_ends[1:]):
                    if region_start <= merged_end:
                        merged_end = max(merged_end, region_end)
                    else:
                        self._collect_insertions_from_region_span(
                            tbx_file,
                            chromosome,
                            merged_start,
                            merged_end,
                            insertion_positions_list,
                            barcode_indices_list,
                        )
                        merged_start = region_start
                        merged_end = region_end

                self._collect_insertions_from_region_span(
                    tbx_file,
                    chromosome,
                    merged_start,
                    merged_end,
                    insertion_positions_list,
                    barcode_indices_list,
                )
                if insertion_positions_list:
                    insertion_positions = np.asarray(insertion_positions_list, dtype=np.int64)
                    barcode_indices = np.asarray(barcode_indices_list, dtype=np.int32)
                    (
                        overlapping_insertion_indices,
                        overlapping_region_indices,
                    ) = region_overlap_lookup.all_overlaps_both(
                        insertion_positions,
                        insertion_positions + 1,
                        np.arange(len(insertion_positions), dtype=np.int64),
                    )
                    np.add.at(
                        chromosome_batch_insertion_counts,
                        (
                            barcode_indices[overlapping_insertion_indices],
                            region_columns[overlapping_region_indices],
                        ),
                        1,
                    )

        return chromosome_batch_insertion_counts

    def _collect_insertions_from_region_span(
        self,
        tbx_file,
        chromosome: str,
        region_span_start: int,
        region_span_end: int,
        insertion_positions_list,
        barcode_indices_list,
    ) -> None:
        """Collect valid fragment endpoints inside one annotation span."""
        # Tabix expects a zero-based query start, while GTF-derived spans are one-based.
        for fragment in tbx_file.fetch(chromosome, max(0, region_span_start - 1), region_span_end):
            fragment_fields = fragment.split("\t", 4)
            barcode_index = self.barcode_to_row.get(fragment_fields[3])
            if barcode_index is None:
                continue

            fragment_start_position = int(fragment_fields[1])
            start_insertion_position = fragment_start_position + 1
            end_insertion_position = int(fragment_fields[2])
            fragment_size = end_insertion_position - fragment_start_position
            if not self.nuc_min_frags <= fragment_size <= self.nuc_max_frags:
                continue

            if region_span_start <= start_insertion_position < region_span_end:
                insertion_positions_list.append(start_insertion_position)
                barcode_indices_list.append(barcode_index)
            if region_span_start <= end_insertion_position < region_span_end:
                insertion_positions_list.append(end_insertion_position)
                barcode_indices_list.append(barcode_index)

    def compute_qc_metrics(
        self,
        nucleosome_df: pd.DataFrame,
        insertion_counts_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Calculate ATAC metrics and return cells that pass the configured filters."""
        start_time = time.time()

        qc_metrics_df = insertion_counts_df.merge(nucleosome_df, on="barcode", how="left")

        nucleosome_cols = [
            "n_frags",
            "n_nucleosome_free_frags",
            "n_mono_frags",
            "n_di_frags",
            "n_multi_frags",
        ]
        qc_metrics_df[nucleosome_cols] = qc_metrics_df[nucleosome_cols].fillna(0).astype(int)

        tss_window_density = qc_metrics_df["window"] / self.window
        flank_density = qc_metrics_df["flank"] / self.norm
        flank_density_normalized = np.maximum(flank_density, self.min_norm)
        qc_metrics_df["tss_enrichment"] = (
            (2 * tss_window_density) / flank_density_normalized
        ).round(3)

        qc_metrics_df["promoter_ratio"] = np.where(
            qc_metrics_df["n_frags"] > 0,
            qc_metrics_df["promoter"] / (qc_metrics_df["n_frags"] * 2),
            np.nan,
        )

        qc_metrics_df["nucleosome_ratio"] = np.where(
            qc_metrics_df["n_nucleosome_free_frags"] > 0,
            (
                qc_metrics_df["n_mono_frags"]
                + qc_metrics_df["n_di_frags"]
                + qc_metrics_df["n_multi_frags"]
            )
            / qc_metrics_df["n_nucleosome_free_frags"],
            np.nan,
        )

        qc_metrics_df["reads_in_tss"] = qc_metrics_df["window"]
        qc_metrics_df["reads_in_promoter"] = qc_metrics_df["promoter"]

        qc_metrics_df = qc_metrics_df.drop(columns=["window", "flank", "promoter"])

        passing_qc_metrics_df = qc_metrics_df[
            (qc_metrics_df.tss_enrichment >= self.min_tss)
            & (qc_metrics_df.n_frags.between(self.min_frags_per_cell, self.max_frags_per_cell))
        ]

        total_cells = len(qc_metrics_df)
        high_quality_cells = len(passing_qc_metrics_df)

        logger.info("Filtering results:")
        logger.info("  Total cells processed: %s", f"{total_cells:,}")
        logger.info("  Cells passing filters: %s", f"{high_quality_cells:,}")

        if high_quality_cells > 0:
            median_tss = passing_qc_metrics_df["tss_enrichment"].median()
            median_frags = passing_qc_metrics_df["n_frags"].median()
            logger.info("  Median TSS (filtered): %.3f", median_tss)
            logger.info("  Median n_frags (filtered): %.0f", median_frags)

            min_passing_frags = passing_qc_metrics_df["n_frags"].min()
            max_passing_frags = passing_qc_metrics_df["n_frags"].max()
            logger.info(
                "  Fragment count range: %s - %s",
                f"{min_passing_frags:,}",
                f"{max_passing_frags:,}",
            )
        else:
            raise ValueError(
                "No cells passed the filtering criteria: "
                f"tss_enrichment >= {self.min_tss}, "
                f"n_frags >= {self.min_frags_per_cell}, "
                f"n_frags <= {self.max_frags_per_cell}"
            )

        elapsed_time = time.time() - start_time
        logger.info(
            "Completed ATAC QC: %s cells passed filters in %.1fs (%.2f min)",
            high_quality_cells,
            elapsed_time,
            elapsed_time / 60,
        )

        passing_qc_metrics_df.columns = [
            "atac_" + col if col != "barcode" else col for col in passing_qc_metrics_df.columns
        ]

        return passing_qc_metrics_df

    def create_tss_regions(
        self,
        transcript_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Create a centered TSS window and two normalization flanks per transcript."""
        regions = list()

        for transcript in transcript_df.drop_duplicates(
            subset=["chrom", "start"], keep="first"
        ).itertuples(index=False):
            tss_start_pos, chromosome = transcript.start, transcript.chrom

            half_width = self.window // 2
            window_start = max(1, tss_start_pos - half_width)
            window_end = window_start + self.window - 1

            upstream_flank_start = max(1, tss_start_pos - self.flank)
            upstream_flank_end = max(1, tss_start_pos - self.flank + self.norm - 1)

            downstream_flank_start = tss_start_pos + self.flank - self.norm + 1
            downstream_flank_end = tss_start_pos + self.flank

            regions.extend(
                [
                    {
                        "chrom": chromosome,
                        "start": window_start,
                        "end": window_end,
                        "type": "window",
                    },
                    {
                        "chrom": chromosome,
                        "start": upstream_flank_start,
                        "end": upstream_flank_end,
                        "type": "flank",
                    },
                    {
                        "chrom": chromosome,
                        "start": downstream_flank_start,
                        "end": downstream_flank_end,
                        "type": "flank",
                    },
                ]
            )

        return pd.DataFrame(regions).sort_values(["chrom", "start"]).reset_index(drop=True)

    def create_promoter_regions(
        self,
        genes_df: pd.DataFrame,
        region_span: tuple[int, int] = (2000, 100),
    ) -> pd.DataFrame:
        """Create strand-aware promoter regions for each gene."""
        upstream, downstream = region_span
        regions = list()

        for gene in genes_df.itertuples(index=False):
            chromosome, strand_orientation = gene.chrom, gene.strand
            tss_start_pos = gene.start if strand_orientation == "+" else gene.end
            if strand_orientation == "+":
                promoter_start = max(1, tss_start_pos - upstream)
                promoter_end = tss_start_pos + downstream
            else:
                promoter_start = max(1, tss_start_pos - downstream)
                promoter_end = tss_start_pos + upstream

            regions.append(
                {
                    "chrom": chromosome,
                    "start": promoter_start,
                    "end": promoter_end,
                    "type": "promoter",
                }
            )

        return pd.DataFrame(regions).sort_values(["chrom", "start"]).reset_index(drop=True)

    def _exclude_scaffold_chromosomes(
        self,
        annotation_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Remove scaffold chromosomes and optionally the mitochondrial chromosome."""
        original_count = len(annotation_df)
        original_chromosomes = set(annotation_df["chrom"])

        annotation_df = annotation_df[
            ~annotation_df["chrom"].str.startswith(
                ("NW_", "NT_", "NG_", "chrGL", "chrUn_", "chrJH", "chrEB")
            )
        ]

        if self.skip_chr_m:
            annotation_df = annotation_df[annotation_df["chrom"] != "chrM"]

        remaining_count = len(annotation_df)
        if original_count != remaining_count:
            remaining_chromosomes = set(annotation_df["chrom"])
            excluded_chromosomes = sorted(original_chromosomes - remaining_chromosomes)

            logger.info(
                "    Excluded chromosomes: %s -> %s (%s removed)",
                original_count,
                remaining_count,
                original_count - remaining_count,
            )
            logger.info("    Excluded chromosomes: %s", excluded_chromosomes)
        return annotation_df

    def _split_gtf_file(
        self,
        gtf_path: str,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Read gene and transcript coordinates from a GTF file."""
        logger.info("Processing Annotation GTF file: %s", gtf_path)

        gtf_column_names = [
            "chrom",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]
        gtf_annotation_df = pd.read_csv(
            gtf_path, sep="\t", comment="#", header=None, names=gtf_column_names
        )

        gtf_annotation_df["chrom"] = gtf_annotation_df["chrom"].apply(
            lambda chromosome: (
                f"chr{chromosome}"
                if not chromosome.startswith("chr")
                and not chromosome.startswith(("NC_", "NW_", "NT_", "NG_"))
                else chromosome
            )
        )

        genes_df = gtf_annotation_df[gtf_annotation_df["feature"] == "gene"][
            ["chrom", "start", "end", "strand"]
        ].copy()
        transcripts_df = gtf_annotation_df[gtf_annotation_df["feature"] == "transcript"][
            ["chrom", "start", "end", "strand"]
        ].copy()

        logger.info("  Excluding scaffold chromosomes from gene features:")
        genes_df = self._exclude_scaffold_chromosomes(genes_df)
        logger.info("  Excluding scaffold chromosomes from transcript features:")
        transcript_df = self._exclude_scaffold_chromosomes(transcripts_df)

        logger.info(
            "    Extracted %s gene records and %s TSS records",
            len(genes_df),
            len(transcript_df),
        )

        return genes_df, transcript_df

    def setup_annotation_regions(self) -> None:
        """Load the GTF file and cache its TSS, flank, and promoter regions."""
        genes_df, transcript_df = self._split_gtf_file(self.annotation_file)

        logger.info("  Creating TSS regions (window=%sbp, flank=%sbp)", self.window, self.flank)
        tss_regions_df = self.create_tss_regions(transcript_df)
        logger.info("    Created %s TSS regions", len(tss_regions_df))

        logger.info("  Creating promoter regions (2000bp upstream, 100bp downstream)")
        promoter_regions_df = self.create_promoter_regions(genes_df, (2000, 100))
        logger.info("    Created %s promoter regions", len(promoter_regions_df))

        self.annotation_regions = (
            pd.concat([tss_regions_df, promoter_regions_df], ignore_index=True)
            .sort_values(["chrom", "start"])
            .reset_index(drop=True)
        )
        self.valid_chromosomes = self.annotation_regions["chrom"].unique()


def run_atac_qc(
    output_path: str,
    annotation_file: str,
    atac_fragments_path: str,
    per_barcode_metrics_path: str,
) -> None:
    """Run the ATAC QC workflow and write `atac_qc.csv`."""
    validate_tabix_index(atac_fragments_path)

    qc_metrics_engine = QCMetricsEngine(
        annotation_file=annotation_file,
        atac_fragments_path=atac_fragments_path,
        per_barcode_metrics_path=per_barcode_metrics_path,
    )

    logger.info("Step 1: Setting up annotation regions")

    qc_metrics_engine.setup_annotation_regions()

    logger.info("Step 2: Classifying nucleosomes and counting fragments")
    nucleosome_df = qc_metrics_engine.nucleosome_classification()

    logger.info("Step 3: Computing TSS enrichment, promoter ratio, and nucleosome ratio metrics")

    qc_metrics_df = qc_metrics_engine.compute_qc_metrics(
        nucleosome_df=nucleosome_df,
        insertion_counts_df=qc_metrics_engine.count_fragments_and_insertions(),
    )
    logger.info("Step 4: Exporting ATAC QC metrics to CSV")

    save_to_csv(
        qc_metrics_df=qc_metrics_df,
        output_path=output_path,
    )


def validate_tabix_index(atac_fragments_path: str) -> None:
    """Require the sibling tabix index needed for random access into fragments."""
    fragment_path = Path(atac_fragments_path)
    index_candidates = list(
        dict.fromkeys(
            [
                fragment_path.with_suffix(fragment_path.suffix + ".tbi"),
                Path(str(fragment_path) + ".tbi"),
            ]
        )
    )
    if not any(index_path.exists() for index_path in index_candidates):
        index_names = ", ".join(str(index_path) for index_path in index_candidates)
        raise FileNotFoundError(
            "ATAC fragments must be bgzip-compressed and tabix-indexed. "
            f"Expected an index at one of: {index_names}"
        )


def save_to_csv(
    qc_metrics_df: pd.DataFrame,
    output_path: str,
) -> None:
    """Write ATAC QC metrics to `atac_qc.csv` in the output directory."""
    output_directory = Path(output_path)
    output_directory.mkdir(parents=True, exist_ok=True)
    output_csv_path = output_directory / "atac_qc.csv"
    qc_metrics_df.to_csv(output_csv_path, index=False)
    logger.info("Results saved to: %s", output_csv_path)


def main() -> None:
    """CLI entry point for ATAC QC."""
    parser = argparse.ArgumentParser(
        description="Compute ATAC QC metrics from portable Cell Ranger ATAC outputs."
    )
    parser.add_argument(
        "--output-path",
        required=True,
        help="Directory where atac_qc.csv will be written.",
    )
    parser.add_argument("--annotation-file", required=True, help="Path to GTF annotation file.")
    parser.add_argument(
        "--atac-fragments-path",
        required=True,
        help="Path to bgzip/tabix-indexed Cell Ranger ATAC fragments.tsv.gz.",
    )
    parser.add_argument(
        "--per-barcode-metrics-path",
        required=True,
        help="Path to Cell Ranger per_barcode_metrics.csv.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()

    run_atac_qc(
        output_path=args.output_path,
        annotation_file=args.annotation_file,
        atac_fragments_path=args.atac_fragments_path,
        per_barcode_metrics_path=args.per_barcode_metrics_path,
    )


if __name__ == "__main__":
    main()
