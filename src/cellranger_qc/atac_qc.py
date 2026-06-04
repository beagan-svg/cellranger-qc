"""
Computes TSS enrichment, promoter ratio, and nucleosome metrics.

IMPORTANT: PLEASE READ FIRST

In this script, an insertion is one endpoint of an ATAC fragment. Biologically,
these endpoints represent Tn5 cut/insertion sites, where the ATAC-seq enzyme cut
accessible DNA and inserted sequencing adapters. Counting insertion points helps
measure how often accessible chromatin occurs near TSS, flank, and promoter
regions.

Coordinate formats:
ATAC fragments: BED-like, 0-based start, exclusive end.
GTF regions: 1-based start, inclusive end.
NCLS: half-open intervals [start, end).

Version 2.0 Update:
For each chromosome, sort its annotation regions, combine overlapping spans for
fewer fragment-file fetches, fetch fragments once per combined span, collect
valid insertion positions, use NCLS to map those insertions back to the exact
original regions, and update the numpy count matrix.

Insertion counting summary:
_count_insertions_for_chromosome_batch coordinates the region-level counting.
For each chromosome, it takes the annotation regions for that chromosome,
combines overlapping or nearby spans so the fragment file is fetched fewer
times, asks _collect_insertions_from_region_span to collect insertion points,
uses region_overlap_lookup to map those insertions back to exact annotation
region types, and adds counts to the barcode-by-region-type matrix.

_collect_insertions_from_region_span collects candidate insertions from the
ATAC fragment file. It fetches fragments overlapping a broader chromosome span,
skips fragments whose barcode is not a valid cell, skips fragments outside the
configured size range, converts each fragment into insertion points
``start + 1`` and ``end``, and keeps only insertion points that fall inside the
span being fetched.

region_overlap_lookup.all_overlaps_both searches for the annotation region that 
each insertion point overlaps.
"""

import argparse
import logging
import os
import time
from concurrent.futures import ProcessPoolExecutor
from typing import Tuple

import numpy as np
import pandas as pd
import polars as pl
import pysam
import warnings
import anndata as ad
from ncls import NCLS

from cellranger_qc import __version__

warnings.filterwarnings("ignore")

MAX_INSERTION_WORKERS = 8

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

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
        """
        Initialize ATAC QC Metrics

        Parameters
        ----------
        annotation_file : str
            Path to the annotation GTF file
        atac_fragments_path : str
            Path to the ATAC fragments TSV file (tab-separated, may be gzipped)
        per_barcode_metrics_path : str
            Path to the per-barcode metrics CSV file
        window : int, optional
            Width of TSS window in base pairs (centered on TSS). Default is 101.
        flank : int, optional
            Distance from TSS to flanking regions in base pairs. Default is 2000.
        norm : int, optional
            Size of normalization region in base pairs for TSS enrichment calculation.
            Default is 100.
        min_norm : float, optional
            Minimum normalization value for flanking region count in TSS enrichment
            calculation. Prevents division by very small numbers. Default is 0.2.
        skip_chr_m : bool, optional
            Whether to skip mitochondrial chromosome (chrM) during processing.
            Default is True.
        min_tss : float, optional
            Minimum TSS enrichment score threshold for filtering cells. Default is 1.0.
        min_frags_per_cell : int, optional
            Minimum number of fragments per cell for filtering. Default is 1000.
        max_frags_per_cell : int, optional
            Maximum number of fragments per cell for filtering. Default is 100000000.
        nuc_min_frags : int, optional
            Minimum fragment size in base pairs for nucleosome classification.
            Default is 10.
        nuc_max_frags : int, optional
            Maximum fragment size in base pairs for nucleosome classification.
            Default is 2000.
        nuc_len : int, optional
            Expected nucleosome length in base pairs used for fragment classification.
            Default is 147.
        """
        self.annotation_file = annotation_file
        self.atac_fragments_path = atac_fragments_path
        per_barcode_metrics_df = pd.read_csv(
            per_barcode_metrics_path,
            usecols=["barcode", "is_cell"],
        )
        
        # Filtered barcodes (is_cell=1)
        self.is_cell_bc = set(per_barcode_metrics_df[per_barcode_metrics_df["is_cell"] == 1]["barcode"])
        total_barcodes = len(per_barcode_metrics_df)
        logging.info(f"Filtered cells: {len(self.is_cell_bc)} out of {total_barcodes} total barcodes "
                    f"({len(self.is_cell_bc)/total_barcodes*100:.1f}%)")
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
        self.barcode_to_row = {
            barcode: row for row, barcode in enumerate(self.barcodes)
        }
        logging.info(
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
        """
        Classify nucleosomes by fragment size for each barcode.
        
        Pseudocode
        ----------
        1. Read the fragment file and keep chromosome, start, end, and barcode.
        2. Keep fragments from valid cell barcodes and valid annotation chromosomes.
        3. Calculate each fragment size as ``end - start``.
        4. Keep fragments within the configured size range.
        5. Convert each fragment size into a nucleosome class.
        6. Count total, mono-, di-, and multi-nucleosome fragments per barcode.

        Data example
        ------------
        With ``nuc_len = 147``:
            barcode  start  end  size  class
            AAAC     100    220  120   mono
            AAAC     300    600  300   di
            TTGG     100    620  520   multi

        The returned counts are:
            barcode  n_frags  n_mono_frags  n_di_frags  n_multi_frags
            AAAC     2        1             1           0
            TTGG     1        0             0           1

        Returns
        -------
        nucleosome_df : pd.DataFrame
            DataFrame with nucleosome counts per barcode. Columns:
                - 'barcode': Cell barcode (str)
                - 'n_frags': Total fragment count (int)
                - 'n_mono_frags': Mono-nucleosome fragments (1x nuc_len) (int)
                - 'n_di_frags': Di-nucleosome fragments (2x nuc_len) (int)
                - 'n_multi_frags': Multi-nucleosome fragments (>=3x nuc_len) (int)
        """
        logging.info(f"Fragments file path: {self.atac_fragments_path}")

        file_size_bytes = os.path.getsize(self.atac_fragments_path)
        file_size_mb = file_size_bytes / (1024 * 1024)
        file_size_gb = file_size_bytes / (1024 * 1024 * 1024)
        if file_size_gb >= 1:
            logging.info(f"ATAC fragment file size: {file_size_gb:.2f} GB ({file_size_bytes:,} bytes)")
        else:
            logging.info(f"ATAC fragment file size: {file_size_mb:.2f} MB ({file_size_bytes:,} bytes)")
        start_time = time.time()

        filtered_fragments_df = (
            pl.scan_csv(
                self.atac_fragments_path,
                has_header=False,
                separator="\t",
                comment_prefix="#",
            )
            # Select only the first 4 columns and assign canonical names/types.
            # Extra columns (if any) are ignored; missing values become nulls and are filtered out.
            .select([
                pl.col("column_1").cast(pl.Utf8).alias("chrom"),
                pl.col("column_2").cast(pl.Int64).alias("start"),
                pl.col("column_3").cast(pl.Int64).alias("end"),
                pl.col("column_4").cast(pl.Utf8).alias("barcode"),
            ])
            .drop_nulls(["chrom", "start", "end", "barcode"])  
            .filter(pl.col("barcode").is_in(self.is_cell_bc))
            .filter(pl.col("chrom").is_in(self.valid_chromosomes))
            .with_columns(
                [
                    (pl.col("end") - pl.col("start")).alias("frag_size"),
                ]
            )
            .filter((pl.col("frag_size") >= self.nuc_min_frags) & (pl.col("frag_size") <= self.nuc_max_frags))
            .with_columns(
                (
                    (pl.col("frag_size") / self.nuc_len)
                    .floor()
                    .cast(pl.Int32)
                    .add(1)
                    .alias("nuc_type")
                )
            )
            .collect()
        )

        nucleosome_df = (
            filtered_fragments_df.group_by(["barcode"])
            .agg([
                pl.len().alias("n_frags"),
                pl.col("nuc_type").filter(pl.col("nuc_type") == 1).count().alias("n_mono_frags"),
                pl.col("nuc_type").filter(pl.col("nuc_type") == 2).count().alias("n_di_frags"),
                pl.col("nuc_type").filter(pl.col("nuc_type") >= 3).count().alias("n_multi_frags"),
            ])
            .to_pandas()
        )
        
        logging.info(f"Total fragments: {filtered_fragments_df.height:,}")
        nucleosome_time = time.time() - start_time
        logging.info(f"Nucleosome counting processing time: {nucleosome_time:.2f} sec ({nucleosome_time/60:.2f} min)")

        return nucleosome_df

    def count_fragments_and_insertions(self) -> pd.DataFrame:
        """
        Count fragment insertions in TSS window, flank, and promoter regions.

        Fragment insertions are counted when fragment start or end positions fall
        within the specified regions. Chromosomes are split across workers so
        each worker gets roughly the same number of annotation regions to process.

        Here, an insertion is one endpoint of an ATAC fragment. Biologically,
        these endpoints represent Tn5 cut/insertion sites in accessible DNA.

        Input flow
        ----------
        GTF annotation file:
            -> create TSS window regions
            -> create flank regions
            -> create promoter regions

        ATAC fragments file:
            -> read fragments overlapping those regions
            -> convert each fragment into insertion points
            -> count which insertion points land inside each region type

        Pseudocode
        ----------
        1. Open the fragment index and keep only annotation chromosomes that exist.
        2. Split chromosomes across workers so each worker receives a similar
           number of annotation regions.
        3. Count insertion overlaps for each chromosome batch in parallel.
        4. Sum the worker count matrices into one barcode-by-region matrix.
        5. Return one row per barcode with nonzero window, flank, or promoter counts.

        Data example
        ------------
        Annotation regions:
            chr1  100  200  window
            chr1  500  600  promoter

        Fragments:
            chr1  120  180  AAAC
            chr1  550  700  TTGG

        The first fragment contributes insertions 121 and 180 to ``window``.
        The second contributes insertion 551 to ``promoter``; 700 is outside
        the promoter. The output counts are:
            barcode  window  flank  promoter
            AAAC     2       0      0
            TTGG     0       0      1

        Returns
        -------
        insertion_counts_df : pd.DataFrame
            DataFrame with insertion counts per barcode and region type. 
            Columns:
                - 'barcode': Cell barcode (str)
                - 'window': Insertions in TSS window (centered on TSS) (int)
                - 'flank': Insertions in TSS flanking regions (int)
                - 'promoter': Insertions in promoter regions (int)
        """
        start_time = time.time()
        insertion_counts = np.zeros((len(self.barcodes), 3), dtype=np.int64)

        with pysam.TabixFile(self.atac_fragments_path) as tbx_file:
            available_contigs = set(tbx_file.contigs)

        chromosome_batches = list()
        for worker_number in range(MAX_INSERTION_WORKERS):
            chromosome_batches.append(list())
        batch_sizes = [0] * MAX_INSERTION_WORKERS
        for chromosome, chromosome_regions in self.annotation_regions.groupby("chrom", sort=False):
            if chromosome not in available_contigs:
                continue

            worker_index = batch_sizes.index(min(batch_sizes))
            chromosome_batches[worker_index].append((chromosome, chromosome_regions))
            batch_sizes[worker_index] += len(chromosome_regions)

        with ProcessPoolExecutor(max_workers=MAX_INSERTION_WORKERS) as executor:
            for chromosome_batch_counts in executor.map(self._count_insertions_for_chromosome_batch, chromosome_batches):
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
        logging.info("TSS and promoter insertion counting summary:")
        logging.info(f"  Window insertions: {insertion_counts[:, 0].sum():,}")
        logging.info(f"  Flank insertions: {insertion_counts[:, 1].sum():,}")
        logging.info(f"  Promoter insertions (start + end points): {insertion_counts[:, 2].sum():,}")
        logging.info(f"  Processing time: {elapsed_time:.1f}s ({elapsed_time/60:.1f} min)")
        return insertion_counts_df

    def _count_insertions_for_chromosome_batch(self, chromosome_batch) -> np.ndarray:
        """
        Count region-overlapping insertion points for one batch of chromosomes.

        Pseudocode
        ----------
        1. Take the regions assigned to this worker, grouped by chromosome.
        2. For each chromosome, read only the fragment-file ranges that cover
           those regions.
        3. For each valid cell fragment, keep its start and end insertion positions.
        4. Check whether each insertion position lands in a TSS window, flank,
           or promoter region.
        5. Add one count for the matching barcode and region type.

        Data example
        ------------
        ``chromosome_batch`` is a list of chromosome-specific region tables:
            [
                (
                    "chr1",
                    chrom  start  end  type
                    chr1   100    200  window
                    chr1   180    260  promoter
                    chr1   500    600  flank
                ),
                (
                    "chr2",
                    chrom  start  end   type
                    chr2   1000   1100  window
                    chr2   1500   1600  promoter
                )
            ]

        Assume the output rows and columns mean:
            row 0 -> barcode AAAC
            row 1 -> barcode TTGG
            column 0 -> window
            column 1 -> flank
            column 2 -> promoter

        If AAAC has two insertions in a TSS window and TTGG has one insertion
        in a promoter, the returned numpy matrix is:

            [[2, 0, 0],
             [0, 0, 1]]

        Parameters
        ----------
        chromosome_batch : list
            List of ``(chromosome, chromosome_regions)`` pairs. Each
            ``chromosome_regions`` DataFrame contains ``start``, ``end``, and
            ``type`` columns for TSS windows, flanks, or promoters.

        Returns
        -------
        np.ndarray
            Barcode-by-region count matrix with columns ordered as window,
            flank, and promoter.
        """
        region_type_to_column = {"window": 0, "flank": 1, "promoter": 2}
        chromosome_batch_insertion_counts = np.zeros((len(self.barcodes), 3), dtype=np.int64)

        with pysam.TabixFile(self.atac_fragments_path) as tbx_file:
            for chromosome, chromosome_regions in chromosome_batch:
                chromosome_regions = chromosome_regions.sort_values(["start", "end"])
                region_starts = chromosome_regions["start"].to_numpy(dtype=np.int64)
                region_ends = chromosome_regions["end"].to_numpy(dtype=np.int64) + 1 # NCLS uses half-open intervals:
                region_rows = np.arange(len(chromosome_regions), dtype=np.int64)
                region_columns = np.asarray(
                    [region_type_to_column[value] for value in chromosome_regions["type"]],
                    dtype=np.int8,
                )
                region_overlap_lookup = NCLS(region_starts, region_ends, region_rows)
                # ATAC fragment start/end insertion positions collected from the fragments file.
                insertion_positions_list = list()
                # Row indexes into self.barcodes, which comes from per_barcode_metrics.csv.
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
                    # Example overlap lookup:
                    #   region 0 = chr1 100-200 window
                    #   region 1 = chr1 180-260 promoter
                    #   insertion_positions = [121, 190, 251]
                    # NCLS treats each insertion as a one-base interval:
                    #   insertion index 0 -> [121, 122)
                    #   insertion index 1 -> [190, 191)
                    #   insertion index 2 -> [251, 252)
                    # Overlaps:
                    #   121 is in region 0
                    #   190 is in region 0 and region 1
                    #   251 is in region 1
                    # Returned indices:
                    #   overlapping_insertion_indices = [0, 1, 1, 2]
                    #   overlapping_region_indices = [0, 0, 1, 1]
                    overlapping_insertion_indices, overlapping_region_indices = region_overlap_lookup.all_overlaps_both(
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
        """
        Collect valid insertion positions inside one region span.

        Pseudocode
        ----------
        1. Fetch fragments from the tabix file over the region span.
        2. Skip fragments whose barcode is not a cell barcode.
        3. Skip fragments outside the configured fragment-size range.
        4. Convert the fragment start to an insertion position by adding 1.
        5. Append start/end insertions that fall inside the region span.

        Data example
        ------------
        For region span ``chr1:100-200`` and fragment ``chr1 120 180 AAAC``:
            fragment start insertion = 121
            fragment end insertion = 180

        Both positions are appended to ``insertion_positions_list`` and the row
        for barcode ``AAAC`` is appended twice to ``barcode_indices_list``.

        Parameters
        ----------
        tbx_file : pysam.TabixFile
            Open tabix handle for the ATAC fragments file.
        chromosome : str
            Chromosome to query.
        region_span_start : int
            Start of the region span.
        region_span_end : int
            End of the region span.
        insertion_positions_list : list
            Output list populated with insertion positions.
        barcode_indices_list : list
            Output list populated with barcode indices matching ``insertion_positions_list``.
        """
        # tbx fetch expects a 0-based start, but region spans are built from GTF, which are 1-based.
        # When we query tabix by region span, tabix returns fragments that overlap the region span, 
        # not fragments whose insertion points are both inside the span.
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
        load_name: str, 
        nucleosome_df: pd.DataFrame,
        insertion_counts_df: pd.DataFrame,
        ) -> pd.DataFrame:
        """
        Calculates TSS enrichment score, promoter ratio, and nucleosome ratio for each
        barcode, then applies filtering criteria based on TSS enrichment and fragment counts.

        Parameters
        ----------
        load_name : str
            Load Name.
        nucleosome_df : pd.DataFrame
            DataFrame from nucleosome_classification() with nucleosome counts per barcode
        insertion_counts_df : pd.DataFrame
            DataFrame from count_fragments_and_insertions() with insertion counts per barcode

        Returns
        -------
        pd.DataFrame
            DataFrame with passing QC metrics, containing columns:
                - 'barcode': Prefixed barcode (str)
                - 'tss_enrichment': TSS enrichment score (float)
                - 'reads_in_tss': Insertions in TSS window (int)
                - 'n_frags': Total fragments (int)
                - 'n_mono_frags': Mono-nucleosome fragments (int)
                - 'n_di_frags': Di-nucleosome fragments (int)
                - 'n_multi_frags': Multi-nucleosome fragments (int)
                - 'nucleosome_ratio': (di + multi) / mono ratio (float)
                - 'reads_in_promoter': Promoter insertions (int)
                - 'promoter_ratio': Promoter insertions / (n_frags * 2) (float)

            Results include only cells meeting:
                - tss_enrichment >= min_tss
                - n_frags >= min_frags_per_cell
                - n_frags <= max_frags_per_cell
        """
        start_time = time.time()

        # Merge insertion counts with nucleosome signal calculations
        qc_metrics_df = insertion_counts_df.merge(nucleosome_df, on='barcode', how='left')
        
        # Fill missing nucleosome signal counts with zeros (for barcodes without fragments)
        nucleosome_cols = ['n_frags', 'n_mono_frags', 'n_di_frags', 'n_multi_frags']
        qc_metrics_df[nucleosome_cols] = qc_metrics_df[nucleosome_cols].fillna(0).astype(int)
        
        # TSS enrichment calculation
        tss_window_density = qc_metrics_df['window'] / self.window
        flank_density = qc_metrics_df['flank'] / self.norm
        flank_density_normalized = np.maximum(flank_density, self.min_norm)
        qc_metrics_df['tss_enrichment'] = ((2 * tss_window_density) / flank_density_normalized).round(3)
        
        # Promoter ratio calculation
        qc_metrics_df['promoter_ratio'] = np.where(
            qc_metrics_df['n_frags'] > 0,
            qc_metrics_df['promoter'] / (qc_metrics_df['n_frags'] * 2),
            np.nan
        )
        
        # Nucleosome ratio calculation
        qc_metrics_df['nucleosome_ratio'] = np.where(
            qc_metrics_df['n_mono_frags'] > 0,
            (qc_metrics_df['n_di_frags'] + qc_metrics_df['n_multi_frags']) / qc_metrics_df['n_mono_frags'],
            np.nan
        )
        
        # Create final column names before dropping temporary columns
        qc_metrics_df['reads_in_tss'] = qc_metrics_df['window']
        qc_metrics_df['reads_in_promoter'] = qc_metrics_df['promoter']
        
        # Drop temporary columns used for calculations
        qc_metrics_df = qc_metrics_df.drop(columns=['window', 'flank', 'promoter'])
        
        passing_qc_metrics_df = qc_metrics_df[
            (qc_metrics_df.tss_enrichment >= self.min_tss) &
            (qc_metrics_df.n_frags.between(
                self.min_frags_per_cell,
                self.max_frags_per_cell
            ))
        ]

        # Log Summary Statistics
        total_cells = len(qc_metrics_df)
        high_quality_cells = len(passing_qc_metrics_df)
        
        logging.info(f"Filtering results:")
        logging.info(f"  Total cells processed: {total_cells:,}")
        logging.info(f"  Cells passing filters: {high_quality_cells:,}")

        if high_quality_cells > 0:
            median_tss = passing_qc_metrics_df['tss_enrichment'].median()
            median_frags = passing_qc_metrics_df['n_frags'].median()
            logging.info(f"  Median TSS (filtered): {median_tss:.3f}")
            logging.info(f"  Median n_frags (filtered): {median_frags:.0f}")
            
            min_passing_frags = passing_qc_metrics_df['n_frags'].min()
            max_passing_frags = passing_qc_metrics_df['n_frags'].max()
            logging.info(f"  Fragment count range: {min_passing_frags:,} - {max_passing_frags:,}")
        else:
            raise ValueError(f"  No cells passed the filtering criteria: tss_enrichment >= {self.min_tss}, n_frags >= {self.min_frags_per_cell}, n_frags <= {self.max_frags_per_cell}")

        elapsed_time = time.time() - start_time
        logging.info(f"Completed {load_name}: {high_quality_cells} cells passed filters in {elapsed_time:.1f}s ({elapsed_time/60:.2f} min)")
        
        # Add "atac_" prefix to all columns except "barcode"
        passing_qc_metrics_df.columns = ['atac_' + col if col != 'barcode' else col for col in passing_qc_metrics_df.columns]
        
        return passing_qc_metrics_df

    def create_tss_regions(
        self, 
        transcript_df: pd.DataFrame,
        ) -> pd.DataFrame:
        """
        Generate TSS window and flanking regions from transcript annotations.

        Creates three region types per TSS:
        1. Window: Centered on TSS with width = self.window (default 101bp)
        2. Upstream flank: 100bp region at TSS - flank to TSS - flank + norm
        3. Downstream flank: 100bp region at TSS + flank - norm to TSS + flank

        Data example
        ------------
        For transcript ``chr1 start=5000 end=9000`` with defaults
        ``window=101``, ``flank=2000``, and ``norm=100``:
            TSS position = 5000
            TSS window = chr1 4950 5050 window
            upstream flank = chr1 3000 3099 flank
            downstream flank = chr1 6901 7000 flank

        Parameters
        ----------
        transcript_df : pd.DataFrame
            DataFrame with columns:
                - 'chrom': Chromosome name (str)
                - 'start': TSS position (int, 1-based)
                - 'end': Gene end position (int, 1-based)
        Returns
        -------
        pd.DataFrame
            DataFrame with columns:
                - 'chrom': Chromosome name (str)
                - 'start': Region start position (int, 1-based)
                - 'end': Region end position (int, 1-based)
                - 'type': Region type ('window' or 'flank') (str)
            Sorted by chrom, then start. Duplicates removed.
        """
        regions = list()
       
        for transcript in transcript_df.drop_duplicates(subset=["chrom", "start"], keep='first').itertuples(index=False):
            tss_start_pos, chromosome = transcript.start, transcript.chrom

            # Create TSS window (101bp centered on TSS)
            half_width = self.window // 2
            window_start = max(1, tss_start_pos - half_width)
            window_end = window_start + self.window - 1

            # Create upstream flank (100bp at TSS - 2000 to TSS - 1901)
            upstream_flank_start = max(1, tss_start_pos - self.flank)
            upstream_flank_end = max(1, tss_start_pos - self.flank + self.norm - 1)

            # Create downstream flank (100bp at TSS + 1901 to TSS + 2000)
            downstream_flank_start = tss_start_pos + self.flank - self.norm + 1
            downstream_flank_end = tss_start_pos + self.flank

            regions.extend([
                {"chrom": chromosome, "start": window_start, "end": window_end, "type": "window"},
                {"chrom": chromosome, "start": upstream_flank_start, "end": upstream_flank_end, "type": "flank"},
                {"chrom": chromosome, "start": downstream_flank_start, "end": downstream_flank_end, "type": "flank"},
            ])

        return pd.DataFrame(regions).sort_values(['chrom', 'start']).reset_index(drop=True)

    def create_promoter_regions(
        self, 
        genes_df: pd.DataFrame, 
        region_span: Tuple[int, int] = (2000, 100)
        ) -> pd.DataFrame:
        """
        Create strand-aware promoter regions from annotations file.

        For each gene, defines promoter region relative to TSS based on strand:
        - Forward strand (+): upstream bp upstream, downstream bp downstream
        - Reverse strand (-): upstream bp downstream, downstream bp upstream

        Data example
        ------------
        With the default ``region_span=(2000, 100)``:

        For ``chr1 start=5000 end=9000 strand=+``:
            TSS position = 5000
            promoter = chr1 3000 5100 promoter

        For ``chr1 start=5000 end=9000 strand=-``:
            TSS position = 9000
            promoter = chr1 8900 11000 promoter

        Parameters
        ----------
        genes_df : pd.DataFrame
            DataFrame with columns:
                - 'chrom': Chromosome name (str)
                - 'start': Gene start position (int, 1-based)
                - 'end': Gene end position (int, 1-based)
                - 'strand': Strand orientation ('+' or '-') (str)
                - Optional: 'gene_id', 'symbol' (not included in output)
        region_span : Tuple[int, int], optional
            (upstream, downstream) tuple defining promoter span in bp.
            Default is (2000, 100) for 2000bp upstream and 100bp downstream.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns:
                - 'chrom': Chromosome name (str)
                - 'start': Promoter start position (int, 1-based)
                - 'end': Promoter end position (int, 1-based)
                - 'type': Always 'promoter' (str)
            Sorted by chrom, then start. Duplicates removed.
        """
        upstream, downstream = region_span
        regions = list()

        for gene in genes_df.itertuples(index=False):
            chromosome, strand_orientation = gene.chrom, gene.strand
            tss_start_pos = gene.start if strand_orientation == "+" else gene.end
            if strand_orientation == "+":
                promoter_start, promoter_end = max(1, tss_start_pos - upstream), tss_start_pos + downstream
            else:
                promoter_start, promoter_end = max(1, tss_start_pos - downstream), tss_start_pos + upstream

            regions.extend([
                {"chrom": chromosome, "start": promoter_start, "end": promoter_end, "type": "promoter"}
            ])

        return pd.DataFrame(regions).sort_values(['chrom', 'start']).reset_index(drop=True)

    def _exclude_scaffold_chromosomes(
        self, annotation_df: pd.DataFrame, 
        ) -> pd.DataFrame:
        """
        Exclude scaffold chromosomes and mitochondrial chromosome (if skip_chr_m is True).
        
        - Scaffold chromosomes: NW_*, NT_*, NG_*, chrGL*, chrUn_*, chrJH*, chrEB*
        - Mitochondrial chromosome (chrM) 

        Parameters
        ----------
        annotation_df : pd.DataFrame
            DataFrame containing chromosome annotations.

        Returns
        -------
        pd.DataFrame
            DataFrame with scaffolds and optionally chrM excluded.
        """
        original_count = len(annotation_df)
        original_chromosomes = set(annotation_df['chrom'])
        
        # Exclude scaffolds
        annotation_df = annotation_df[~annotation_df['chrom'].str.startswith(('NW_', 'NT_', 'NG_', 'chrGL', 'chrUn_', 'chrJH', 'chrEB'))]
        
        # Exclude chrM if enabled
        if self.skip_chr_m:
            annotation_df = annotation_df[annotation_df['chrom'] != 'chrM']
        
        remaining_count = len(annotation_df)
        if original_count != remaining_count:
            remaining_chromosomes = set(annotation_df['chrom'])
            excluded_chromosomes = sorted(original_chromosomes - remaining_chromosomes)
            
            logging.info(f"    Excluded chromosomes: {original_count} -> {remaining_count} "
                        f"({original_count - remaining_count} removed)")
            logging.info(f"    Excluded chromosomes: {excluded_chromosomes}")
        return annotation_df

    def _split_gtf_file(
        self, 
        gtf_path: str
        ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Process GTF file and organize into gene and transcript annotations.

        Reads GTF file, filters for 'gene' and 'transcript' features, standardizes
        chromosome naming, excludes scaffolds, and returns separate DataFrames.

        Parameters
        ----------
        gtf_path : str
            Path to GTF annotation file (may be gzipped). File should be tab-separated
            with standard GTF columns: chrom, source, feature, start, end, score,
            strand, frame, attribute.

        Returns
        -------
        gene_df : pd.DataFrame
            DataFrame with columns ['chrom', 'start', 'end', 'strand'] for gene features.
            Scaffolds excluded, chromosomes standardized.
        transcript_df : pd.DataFrame
            DataFrame with columns ['chrom', 'start', 'end', 'strand'] for transcript features.
            Scaffolds excluded, chromosomes standardized. Start position represents TSS.
        """
        logging.info(f"Processing Annotation GTF file: {gtf_path}")

        gtf_column_names = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        gtf_annotation_df = pd.read_csv(gtf_path, sep='\t', comment='#', header=None, names=gtf_column_names)
        
        # Standardize chromosome naming
        gtf_annotation_df['chrom'] = gtf_annotation_df['chrom'].apply(
            lambda chromosome: f'chr{chromosome}' if not chromosome.startswith('chr') and not chromosome.startswith(('NC_', 'NW_', 'NT_', 'NG_')) else chromosome
        )
        
        # Split gtf by features into gene and transcript dataframes
        genes_df = gtf_annotation_df[gtf_annotation_df['feature'] == 'gene'][['chrom', 'start', 'end', 'strand']].copy()
        transcripts_df = gtf_annotation_df[gtf_annotation_df['feature'] == 'transcript'][['chrom', 'start', 'end', 'strand']].copy()

        logging.info("  Excluding scaffold chromosomes from gene features:")
        genes_df = self._exclude_scaffold_chromosomes(genes_df)
        logging.info("  Excluding scaffold chromosomes from transcript features:")
        transcript_df = self._exclude_scaffold_chromosomes(transcripts_df)

        logging.info(f"    Extracted {len(genes_df)} gene records and {len(transcript_df)} TSS records")

        return genes_df, transcript_df

    def setup_annotation_regions(self) -> None:
        """
        Load annotation file and create TSS/promoter regions.

        Processes GTF annotation file, extracts gene and transcript features, creates
        TSS window/flank regions and promoter regions.

        Method Effects
        ------------
        Sets self.annotation_regions and self.valid_chromosomes.
        """

        genes_df, transcript_df = self._split_gtf_file(self.annotation_file)

        logging.info(f"  Creating TSS regions (window={self.window}bp, flank={self.flank}bp)")
        tss_regions_df = self.create_tss_regions(transcript_df)
        logging.info(f"    Created {len(tss_regions_df)} TSS regions")

        logging.info(f"  Creating promoter regions (2000bp upstream, 100bp downstream)")
        promoter_regions_df = self.create_promoter_regions(genes_df, (2000, 100))
        logging.info(f"    Created {len(promoter_regions_df)} promoter regions")

        self.annotation_regions = (
            pd.concat([tss_regions_df, promoter_regions_df], ignore_index=True)
            .sort_values(["chrom", "start"])
            .reset_index(drop=True)
        )
        self.valid_chromosomes = self.annotation_regions["chrom"].unique()

def run_atac_qc(
    output_path: str,
    load_name: str,
    organism_common_name: str,
    batch_name_from_vendor: str,
    library_prep_name: str,
    alignment_fs_id: str,
    annotation_file: str,
    atac_fragments_path: str,
    per_barcode_metrics_path: str
    ):
    """
    Run ATAC QC, merge metrics into GEX AnnData, and save a CSV.

    Parameters
    ----------
    output_path : str
        Path to the output directory for saving results
    load_name : str
        Load name
    organism_common_name : str
        Organism common name from the pipeline interface.
    batch_name_from_vendor : str
        Batch name from vendor from the pipeline interface.
    library_prep_name : str
        Library prep name
    alignment_fs_id : str
        Alignment file store ID
    annotation_file : str
        Path to the annotation GTF file
    atac_fragments_path : str
        Path to the ATAC fragments TSV (possibly gzipped)
    per_barcode_metrics_path : str
        Path to the per-barcode metrics CSV
    """

    qc_metrics_engine = QCMetricsEngine(
        annotation_file=annotation_file,
        atac_fragments_path=atac_fragments_path,
        per_barcode_metrics_path=per_barcode_metrics_path
    )
    
    logging.info("=" * 60)
    logging.info("Step 1: Setting up annotation regions")
    logging.info("=" * 60)

    qc_metrics_engine.setup_annotation_regions()

    logging.info("=" * 60)
    logging.info("Step 2: Classifying nucleosomes and counting fragments")
    logging.info("=" * 60)
    nucleosome_df = qc_metrics_engine.nucleosome_classification()

    logging.info("=" * 60)
    logging.info("Step 3: Computing TSS enrichment, promoter ratio, and nucleosome ratio metrics")
    logging.info("=" * 60)

    qc_metrics_df = qc_metrics_engine.compute_qc_metrics(
        load_name=load_name,
        nucleosome_df=nucleosome_df,
        insertion_counts_df=qc_metrics_engine.count_fragments_and_insertions(),
    )
     
    logging.info("=" * 60)
    logging.info("Step 4: Merging ATAC QC metrics into gex AnnData object and exporting to CSV")
    logging.info("=" * 60)
    
    read_and_merge_gex_h5ad(
        qc_metrics_df=qc_metrics_df,
        library_prep_name=library_prep_name,
        alignment_fs_id=alignment_fs_id,
        output_path=output_path
    )
    
    qc_metrics_df = add_metadata(
        qc_metrics_df=qc_metrics_df,
        output_path=output_path,
        alignment_fs_id=alignment_fs_id
    )
    
    save_to_csv(
        qc_metrics_df=qc_metrics_df,
        output_path=output_path,
        alignment_fs_id=alignment_fs_id
    )
    
def read_and_merge_gex_h5ad(
    qc_metrics_df: pd.DataFrame,
    library_prep_name: str,
    alignment_fs_id: str,
    output_path: str
    ) -> None:
    """
    Merge ATAC QC metrics into the GEX AnnData obs table.
    
    Parameters
    ----------
    qc_metrics_df : pd.DataFrame
        DataFrame containing ATAC QC results with prefixed barcodes
    library_prep_name : str
        Library prep name used to construct sample_id format
    alignment_fs_id : str
        Alignment file store ID used to construct sample_id format
    output_path : str
        Path to the output directory containing the gex AnnData file
    """  
     
    qc_metrics_df = qc_metrics_df.set_index(
        qc_metrics_df['barcode'].str.replace(r'-\d+$', '', regex=True)
        + '-' + library_prep_name + '-' + alignment_fs_id
    )
    
    gex_h5ad_path = os.path.join(output_path, f"gex_{alignment_fs_id}.h5ad")
    
    gex_adata = ad.read_h5ad(gex_h5ad_path)
    
    merged_cell_metadata = gex_adata.obs.join(qc_metrics_df, how='left')
     
    unmatched_atac_cells = len(set(qc_metrics_df.index) - set(gex_adata.obs.index))
    if unmatched_atac_cells > 0:
        raise ValueError(f'Found "{unmatched_atac_cells}" ATAC cells not in AnnData ({(unmatched_atac_cells/len(qc_metrics_df)*100):.1f}%)')
        
    # Update the AnnData with ATAC QC metrics
    gex_adata.obs = merged_cell_metadata

    logging.info(f"Saving updated AnnData object to: {gex_h5ad_path}")
    gex_adata.write(gex_h5ad_path)
    logging.info("AnnData object successfully updated with ATAC QC metrics")
    
def add_metadata(
    qc_metrics_df: pd.DataFrame,
    output_path: str,
    alignment_fs_id: str
    ) -> pd.DataFrame:
    """
    Add the cell_member column from the samp.dat file.
    
    Parameters
    ----------
    qc_metrics_df : pd.DataFrame
        DataFrame containing QC metrics results
    output_path : str
        Path to the output directory containing the samp.dat file
    alignment_fs_id : str
        Alignment file store ID used to construct samp.dat filename
    
    Returns
    -------
    pd.DataFrame
        DataFrame with added cell_member column
    """
    sample_metadata_df = pd.read_csv(f"{output_path}/samp.dat_{alignment_fs_id}.csv")
    if "bc" in sample_metadata_df.columns:
        barcode_column = "bc"
    elif "barcodes" in sample_metadata_df.columns:
        barcode_column = "barcodes"
    else:
        raise ValueError(
            f"samp.dat_{alignment_fs_id}.csv must include either a 'bc' or 'barcodes' column"
        )
    
    qc_metrics_df['raw_barcode'] = qc_metrics_df['barcode'].str.replace(r'-\d+$', '', regex=True)
    
    # Merge with samp.dat to add cell_member column (merge on raw_barcode -> barcodes)
    # Cell Member ID is the match identifier for CODE
    qc_metrics_df = qc_metrics_df.merge(
        sample_metadata_df[[barcode_column, 'cell_member']],
        left_on='raw_barcode',
        right_on=barcode_column,
        how='left',
    )
    
    # Drop the temporary columns used for merging
    qc_metrics_df = qc_metrics_df.drop(columns=['raw_barcode', barcode_column])
    
    return qc_metrics_df

def save_to_csv(
    qc_metrics_df: pd.DataFrame,
    output_path: str,
    alignment_fs_id: str
    ) -> None:
    """
    Save ATAC QC metrics results to CSV file.

    Parameters
    ----------
    qc_metrics_df : pd.DataFrame
        DataFrame containing ATAC QC metrics results
    output_path : str
        Path to the output directory for saving results
    alignment_fs_id : str
        Alignment file store ID used to construct filename
    """
    output_csv_path = os.path.join(output_path, f"atac_qc_{alignment_fs_id}.csv")
    qc_metrics_df.to_csv(output_csv_path, index=False)
    logging.info(f"Results saved to: {output_csv_path}")


def main():
    """CLI entry point for ATAC QC."""
    parser = argparse.ArgumentParser(
        description="Compute ATAC QC metrics from Cell Ranger ARC/multiome outputs."
    )
    parser.add_argument("--output-path", required=True, help="Directory containing gex_<alignment_fs_id>.h5ad and samp.dat_<alignment_fs_id>.csv")
    parser.add_argument("--load-name", required=True, help="Load name used in QC logs.")
    parser.add_argument("--organism-common-name", default="", help="Organism common name; accepted for pipeline compatibility.")
    parser.add_argument("--batch-name-from-vendor", default="", help="Batch name from vendor; accepted for pipeline compatibility.")
    parser.add_argument("--library-prep-name", required=True, help="Library prep name used to construct sample IDs.")
    parser.add_argument("--alignment-fs-id", required=True, help="Alignment file-store ID, such as AR123.")
    parser.add_argument("--annotation-file", required=True, help="Path to GTF annotation file.")
    parser.add_argument("--atac-fragments-path", required=True, help="Path to bgzip/tabix-indexed ATAC fragments TSV.")
    parser.add_argument("--per-barcode-metrics-path", required=True, help="Path to Cell Ranger per_barcode_metrics.csv.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()

    run_atac_qc(
        output_path=args.output_path,
        load_name=args.load_name,
        organism_common_name=args.organism_common_name,
        batch_name_from_vendor=args.batch_name_from_vendor,
        library_prep_name=args.library_prep_name,
        alignment_fs_id=args.alignment_fs_id,
        annotation_file=args.annotation_file,
        atac_fragments_path=args.atac_fragments_path,
        per_barcode_metrics_path=args.per_barcode_metrics_path,
    )


if __name__ == "__main__":
    main()
