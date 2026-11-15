# fatq-vcf

Germline **WGS / WES** pipeline on **GRCh38**, single-sample:

FASTQ → fastp QC → bwa-mem3 alignment → samtools markdup → GATK BQSR → **HaplotypeCaller** (GVCF, DRAGEN mode) → **mitochondrial** calling (GATK Mutect2 mito pipeline) → **SV/CNV** (Manta + CNVpytor for WGS, gCNV + Manta for WES) → MultiQC.

Built with Nextflow DSL2; every tool runs in a version-pinned container. Design rationale, tool choices and the phase plan live in [docs/PIPELINE_PLAN.md](docs/PIPELINE_PLAN.md).

> **Status:** under construction — see the phase table in the plan. Not yet validated; do not use for clinical reporting.

## Requirements

- Linux node; sized for 80 cores / 125 GB RAM / ≥ 1 TB disk (adjust `--max_cpus`, `--max_memory`)
- Docker (user in the `docker` group) or Apptainer
- Nextflow ≥ 25.04 + Java 17 — e.g. `micromamba create -n nf -c conda-forge -c bioconda nextflow nf-test openjdk=17`

## Quick start

```bash
# 0. One-off: build the local alignment image (bwa-mem3 + bwa-mem2 + samtools)
docker build --network=host -t fatqvcf/align:0.1.0 containers/align

# 1. One-off: reference bundle (~15 GB download, ~2 h for the bwa index) + HG002 chr20 test set
nextflow run main.nf -profile docker --prepare_references

# 2. Run
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results

# Regression test (HG002 chr20 + chrM, minutes)
nextflow run main.nf -profile docker,test
```

### Samplesheet

```csv
sample,lane,fastq_1,fastq_2,seq_type,capture_bed,sex
S1,L001,/data/S1_L001_R1.fq.gz,/data/S1_L001_R2.fq.gz,wgs,,female
S2,L001,/data/S2_R1.fq.gz,/data/S2_R2.fq.gz,wes,/kits/twist_exome_2.0.hg38.bed,male
```

Multiple rows with the same `sample` (different `lane`) are merged after alignment. `capture_bed` is required for `wes`.

### Key options

| Option | Default | Meaning |
|---|---|---|
| `--aligner` | `bwamem3` | `bwamem3` (`--compat=bwa-mem2`, byte-identical to bwa-mem2 output) or `bwamem2` |
| `--skip_bqsr` | `false` | skip base-quality recalibration |
| `--skip_mito` / `--skip_sv_cnv` | `false` | skip those branches |
| `--scatter_count` | `20` | interval shards for BQSR + HaplotypeCaller |
| `--markdup_optical_dist` | `2500` | `samtools markdup -d` (2500 patterned flow cells, 100 HiSeq) |
| `--gcnv_model_dir` | `null` | WES panel-of-normals model → enables GATK gCNV for `wes` samples |

Tool versions are pinned in one place: [conf/containers.config](conf/containers.config).

## Outputs (`results/<sample>/`)

`alignment/*.cram` · `qc/` (fastp, mosdepth, metrics, contamination, sex check) · `variants/*.g.vcf.gz`, `*.filtered.vcf.gz` · `mito/` · `sv_cnv/` · `results/multiqc/`

## Node notes

This node's Docker bridge network has no outbound route, so the `docker` profile runs containers with `--network=host`. Only the reference-download steps need network access.
