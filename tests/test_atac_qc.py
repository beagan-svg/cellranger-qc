import pandas as pd
import pytest

from cellranger_qc.atac_qc import QCMetricsEngine, save_to_csv, validate_tabix_index


class TestAtacOutput:
    """Tests for writing ATAC QC output files."""

    def test_save_to_csv_writes_standalone_output_name(self, tmp_path):
        qc_metrics_df = pd.DataFrame({"barcode": ["AAAC-1"], "atac_tss_enrichment": [6.0]})

        save_to_csv(qc_metrics_df, str(tmp_path))

        output_path = tmp_path / "atac_qc.csv"
        assert output_path.exists()
        assert pd.read_csv(output_path).to_dict("records") == [
            {"barcode": "AAAC-1", "atac_tss_enrichment": 6.0}
        ]


class TestTabixValidation:
    """Tests for ATAC fragments tabix index validation."""

    def test_validate_tabix_index_accepts_sibling_tbi(self, tmp_path):
        fragments_path = tmp_path / "fragments.tsv.gz"
        fragments_path.touch()
        (tmp_path / "fragments.tsv.gz.tbi").touch()

        validate_tabix_index(str(fragments_path))

    def test_validate_tabix_index_requires_sibling_tbi(self, tmp_path):
        fragments_path = tmp_path / "fragments.tsv.gz"
        fragments_path.touch()

        with pytest.raises(FileNotFoundError, match="tabix-indexed"):
            validate_tabix_index(str(fragments_path))


class TestComputeQcMetrics:
    """Tests for ATAC QC metric calculations."""

    def test_compute_qc_metrics_uses_nucleosome_free_denominator(self):
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

        qc_metrics_df = engine.compute_qc_metrics(nucleosome_df, insertion_counts_df)

        assert qc_metrics_df.loc[0, "atac_nucleosome_ratio"] == 2
        assert qc_metrics_df.loc[0, "atac_n_nucleosome_free_frags"] == 2

    def test_compute_qc_metrics_raises_when_no_cells_pass_filters(self):
        engine = object.__new__(QCMetricsEngine)
        engine.window = 101
        engine.norm = 100
        engine.min_norm = 0.2
        engine.min_tss = 5
        engine.min_frags_per_cell = 100
        engine.max_frags_per_cell = 1000

        nucleosome_df = pd.DataFrame(
            {
                "barcode": ["AAAC-1"],
                "n_frags": [10],
                "n_nucleosome_free_frags": [2],
                "n_mono_frags": [1],
                "n_di_frags": [0],
                "n_multi_frags": [0],
            }
        )
        insertion_counts_df = pd.DataFrame(
            {
                "barcode": ["AAAC-1"],
                "window": [10],
                "flank": [100],
                "promoter": [2],
            }
        )

        with pytest.raises(ValueError, match="No cells passed"):
            engine.compute_qc_metrics(nucleosome_df, insertion_counts_df)


class TestAnnotationRegions:
    """Tests for creating annotation regions used by ATAC QC."""

    def test_create_tss_regions_deduplicates_and_builds_expected_windows(self):
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

        tss_regions_df = engine.create_tss_regions(transcript_df)

        assert tss_regions_df.to_dict("records") == [
            {"chrom": "chr1", "start": 3000, "end": 3099, "type": "flank"},
            {"chrom": "chr1", "start": 4950, "end": 5050, "type": "window"},
            {"chrom": "chr1", "start": 6901, "end": 7000, "type": "flank"},
        ]

    def test_create_promoter_regions_is_strand_aware(self):
        engine = object.__new__(QCMetricsEngine)
        genes_df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [5000, 5000],
                "end": [9000, 9000],
                "strand": ["+", "-"],
            }
        )

        promoter_regions_df = engine.create_promoter_regions(genes_df, region_span=(2000, 100))

        assert promoter_regions_df.to_dict("records") == [
            {"chrom": "chr1", "start": 3000, "end": 5100, "type": "promoter"},
            {"chrom": "chr1", "start": 8900, "end": 11000, "type": "promoter"},
        ]

    def test_exclude_scaffold_chromosomes_respects_chr_m_setting(self):
        annotation_df = pd.DataFrame({"chrom": ["chr1", "chrM", "chrUn_1", "NW_123"]})
        engine = object.__new__(QCMetricsEngine)

        engine.skip_chr_m = True
        assert engine._exclude_scaffold_chromosomes(annotation_df)["chrom"].tolist() == ["chr1"]

        engine.skip_chr_m = False
        assert engine._exclude_scaffold_chromosomes(annotation_df)["chrom"].tolist() == [
            "chr1",
            "chrM",
        ]
