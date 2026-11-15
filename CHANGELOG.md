# Changelog

All notable changes to fatq-vcf are documented here. Format: [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added
- **Phase 0 — scaffold**: Nextflow DSL2 pipeline in nf-core layout; Docker + Apptainer profiles;
  version-pinned containers (`conf/containers.config`); local `fatqvcf/align` (fastp + bwa-mem3 +
  bwa-mem2 + samtools) and `fatqvcf/mito` (gatk4 + bwa-mem2 + samtools) images.
- **Phase 1 — references**: `--prepare_references` builds GRCh38 no-alt+decoy analysis set,
  `.fai`/`.dict`/bwa-mem2 index, DRAGstr STR table, GATK hg38 bundle (dbSNP bgzipped, Mills,
  known-indels, WGS calling regions), Broad chrM resources (+ indexes), somalier sites, and a
  HG002 chr20+chrM regression subset with GIAB v4.2.1 chr20 truth.
- **Phase 2 — QC + alignment**: fastp (light trim) streamed into bwa-mem3 `--compat=bwa-mem2`
  (bwa-mem2 fallback) → samtools fixmate/sort/markdup → CRAM; mosdepth, samtools stats,
  VerifyBamID2, somalier, Picard HsMetrics (WES); sample-level QC gates (`bin/qc_gates.py`).
- **Phase 3 — germline SNV/indel**: shared interval scatter; BQSR (scatter/gather); DRAGstr model;
  HaplotypeCaller GVCF `--dragen-mode` (scatter) → MergeVcfs → GenotypeGVCFs → VariantFiltration
  (DRAGEN QUAL hard filter) → bcftools stats.
- **Phase 4 — mitochondria**: GATK shifted-reference pipeline (chrM + shifted-by-8000 realignment,
  Mutect2 `--mitochondria-mode`, liftover/merge/filter, blacklist mask) + haplocheck + Haplogrep 3.
- **Phase 5 — SV/CNV**: Manta germline (`--exome` for WES); CNVpytor read-depth CNV (WGS);
  GATK gCNV case mode (WES, gated on a panel of normals); optional AnnotSV.
- **Phase 6 — validation**: `--validate` runs hap.py (small variants) and truvari (SV) against GIAB;
  `ALIGNER_CONCORDANCE` asserts bwa-mem3 (compat) ≡ bwa-mem2 alignment records.
- **Phase 7 — hardening**: `--pon_build` (WES gCNV panel of normals), nf-test suite, param schema
  (`nextflow_schema.json`), CHANGELOG, MIT LICENSE.

### Notes
- bwa-mem3 is pre-1.0 and pinned to 0.17.0; re-run the Phase 6 aligner head-to-head on any bump.
- Not yet validated end-to-end against GIAB; do not use for clinical reporting until Phase 6 passes.
