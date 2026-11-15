# CLAUDE.md — fatq-vcf

Germline WGS/WES pipeline (FASTQ → SNV/indel + mito + CNV/SV VCFs on GRCh38). Nextflow DSL2,
nf-core layout, every tool in a version-pinned container. **Design + rationale: `docs/PIPELINE_PLAN.md`.**

## Repo map
- `main.nf` — entry; `--prepare_references`, `--validate`, `--pon_build` switch the entry workflow.
- `workflows/` — `fatqvcf.nf` (main), `prepare_references.nf`, `validate.nf`, `pon_build.nf`.
- `subworkflows/local/` — `qc_align`, `prepare_intervals`, `bqsr`, `germline_snv`, `mito`, `sv_cnv`.
- `modules/local/` — one process per file, grouped by tool.
- `conf/` — `base` (resource labels), `containers` (**the only place versions are pinned**),
  `references` (bundle paths), `modules` (per-process args + publishDir), `test*`.
- `bin/` — helper scripts on `$PATH` inside jobs (`qc_gates.py`).
- `lib/samplesheet.nf` — samplesheet parsing/validation.

## Conventions
- **Strict parser is on** (`nextflow.enable.strict`). No top-level `def` in scripts; no C-style
  `(int)` casts — use `as int` / `intdiv`; no `-entry`, use the `--<name>` entry switches.
- Lint before running: `nextflow lint main.nf workflows subworkflows modules lib nextflow.config conf`.
- Containers: add a `withName` mapping in `conf/containers.config`; GATK/Picard match `GATK4_.*`.
- Per-sample `meta` map carries `id, sample, lane, seq_type(wgs|wes), capture_bed, sex, single_end`.
- Keep CRAM (not BAM) and bgzipped+tabixed VCFs throughout.

## Environment (this node)
- Run tools via `export PATH=/home/labs/micromamba/envs/nf/bin:$PATH` (nextflow + nf-test + Java 17).
- Docker needs the group shim in non-login shells: wrap docker commands with `sg docker -c "..."`.
- Docker bridge has **no outbound route** → the `docker` profile uses `--network=host` (set already).
- Build local images once:
  `docker build --network=host -t fatqvcf/align:0.1.0 containers/align`
  `docker build --network=host -t fatqvcf/mito:0.1.0 containers/mito`

## Run
```bash
export PATH=/home/labs/micromamba/envs/nf/bin:$PATH
sg docker -c "nextflow run main.nf -profile docker --prepare_references"     # one-off
sg docker -c "nextflow run main.nf -profile docker,test"                      # chr20 regression
sg docker -c "nextflow run main.nf -profile docker --input samplesheet.csv"   # real run
nf-test test tests/lib/samplesheet.nf.test                                    # fast unit test
```

## Hard rules
- **Do not bump bwa-mem3** (pinned 0.17.0) without re-running the Phase 6 aligner head-to-head
  (`ALIGNER_CONCORDANCE` + hap.py). It is pre-1.0 with breaking changes between minor releases.
- Pushing to GitHub needs the user's own credentials (no `gh`, no stored token here).
