# fatq-vcf — Germline Pipeline Design Plan

**Scope:** FASTQ → QC → GRCh38 alignment → duplicate marking → BQSR → HaplotypeCaller (GVCF) → mitochondrial calling → CNV/SV calling, for **WGS and WES**, **single-sample germline**, orchestrated with **Nextflow**.

**Target machine:** 80 cores, 125 GB RAM, 1.1 TB disk, no GPU, Docker 26 installed.

**Design goals, in priority order:** (1) accuracy against GIAB truth sets, (2) wall-clock speed on one 80-core node, (3) reproducibility (pinned containers, resumable runs).

---

## 1. Tool selection

Every tool below is open-source, containerized (BioContainers / Broad images), and current as of this plan. Pin exact versions in `conf/modules.config` at install time — "latest" below means *check the release page once and pin*.

| Stage | Tool | Why this one |
|---|---|---|
| Read QC + trimming | **fastp** 1.x | Single pass: adapter detection, quality trimming, poly-G trimming (NovaSeq/NextSeq), duplication estimate, HTML+JSON report. Replaces FastQC + Trimmomatic/cutadapt and is ~5× faster. Keep trimming *light* (`--trim_poly_g`, adapter auto-detect, no hard quality cut) — BWA soft-clips, and aggressive trimming hurts indel calling. |
| Optional raw QC | FastQC 0.12 | Only if the lab wants the familiar report; off by default. |
| Alignment | **bwa-mem3** (Fulcrum Genomics fork of bwa-mem2; pin ≥ v0.17.0) in **`--compat=bwa-mem2`** mode, with **bwa-mem2 2.3** as `--aligner bwamem2` fallback | bwa-mem3 is a drop-in replacement for bwa-mem2 that fixes ~20 upstream bugs and benchmarks ~5× faster (README: 93 s → 16 s on 5 M WGS pairs, 16 threads; figure has confounds). In `--compat=bwa-mem2` mode its alignment records are byte-identical to bwa-mem2 2.2.1, so every GATK/Manta assumption holds; it reuses the bwa-mem2 index as-is. **Caveats:** pre-1.0, breaking changes between minor releases, no published GIAB variant-level benchmark, not yet in nf-core — hence compat mode, exact version pin, no `--fast`, and a Phase 6 head-to-head against bwa-mem2 before any of those are relaxed. `shm` (shared-memory index) saves the ~11 GB index load per sample on multi-sample runs. Index build ~72–90 GB RAM once (fits), ~11–20 GB disk; denser SA option costs +12 GB RAM at run time. Flags: `-K 100000000` (deterministic), `-Y` (soft-clip supplementary; required by Manta/GATK best practice), `-R` read group. |
| Sort + duplicates | **samtools** 1.22 (`fixmate -m` → `sort` → `markdup`) | Fully streaming from the aligner — no intermediate SAM/BAM on disk. `markdup -d 2500 -S` marks optical duplicates for patterned flow cells and marks supplementary reads. Output CRAM (lossless, ~50 % of BAM). Alternative if GATK-canonical output is required: `GATK MarkDuplicatesSpark` (slower, same result). |
| BQSR | **GATK** 4.6.x `BaseRecalibrator` / `ApplyBQSR` | Scattered across interval shards and gathered (`GatherBQSRReports`). Known sites: dbSNP, Mills+1000G gold indels, known_indels (GATK hg38 bundle). WES: restricted to capture BED ±100 bp. Exposed as `--skip_bqsr` (see §6). |
| STR calibration | GATK `ComposeSTRTableFile` (once) + `CalibrateDragstrModel` (per sample) | Feeds HaplotypeCaller's **DRAGEN mode** — the main recent accuracy gain in GATK for indels in short tandem repeats. |
| SNV/indel calling | **GATK HaplotypeCaller** `-ERC GVCF --dragen-mode true` | GVCF keeps the door open to joint genotyping later. Scattered over ~20 shards × 4 pairHMM threads. Annotations `-G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation`. |
| Genotyping + filtering | GATK `GenotypeGVCFs` → `VariantFiltration` (DRAGEN hard filters: `QUAL < 10.4139`) | Single sample ⇒ no VQSR. DRAGEN-mode hard filter is the Broad-recommended replacement. Optional second filter path: `NVScoreVariants` + `FilterVariantTranches` (the CNN successor) for a tranche-based score. |
| Mitochondria | **GATK Mitochondria pipeline** (Mutect2 `--mitochondria-mode`, shifted-reference realignment, `FilterMutectCalls --mitochondria-mode`) + **haplocheck** (contamination) + **Haplogrep 3** (haplogroup) | Broad's production mtDNA workflow; handles the circular genome break, heteroplasmy down to ~1–5 % VAF, NuMT filtering via autosomal coverage. |
| SV / CNV — WGS | **Manta** 1.6 (split/paired-read SVs ≥ 50 bp) + **CNVpytor** 1.3 (read-depth CNVs, kb–Mb, with BAF support) | Both work on a single sample with **no panel of normals**. Complementary: Manta gives base-pair-precise breakpoints; CNVpytor catches large/low-mappability events. Merge/annotate with **AnnotSV** (optional). |
| CNV — WES | **GATK GermlineCNVCaller** (case mode vs lab **panel of normals**) + Manta `--exome` | Exome read-depth CNV calling is unreliable without a reference set. PoN is built once per capture kit from ≥ 20–30 lab samples (separate `pon` workflow). Until the PoN exists, gCNV is skipped and Manta-only results are flagged *low-sensitivity*. |
| Alignment / coverage QC | **mosdepth** 0.3 (coverage, fast), `samtools stats/flagstat`, GATK/Picard `CollectWgsMetrics` or `CollectHsMetrics` (WES), **VerifyBamID2** (contamination), **somalier** (sex check, sample swap) | Clinically essential gates: contamination > 2–3 % and sex mismatch should fail a sample. |
| Report | **MultiQC** 1.3x | One HTML per run aggregating fastp, samtools, mosdepth, GATK metrics, VerifyBamID, bcftools stats. |
| Orchestration | **Nextflow** 25.x, DSL2, nf-core module layout, **nf-test** | Resume, per-process containers, native scatter/gather, `-with-report/-with-trace` for tuning. |

**Considered and not chosen**
- *NVIDIA Parabricks* — fastest option (~1 h/WGS) but needs a GPU. Revisit if one is added; the Nextflow process boundaries make swapping `fp_bwa`/`haplotypecaller` in trivial.
- *DRAGMAP* (Illumina's open-source DRAGEN mapper) — pairs with `--dragen-mode`, but upstream is stale since 2022; bwa-mem2/3 are better maintained and the STR model is aligner-agnostic.
- *minibwa* (June 2026 preprint, ~2× bwa-mem2) — newer still, less tested than bwa-mem3, and bwa-mem3's own benchmark beats it. Watch, don't adopt.
- *bwa-mem3 `--fast`* — a further ~1.7× but not alignment-identical ("trades some sensitivity/specificity at the extremes"). Off until Phase 6 shows no GIAB F1 loss.
- *DeepVariant* — typically beats HaplotypeCaller on GIAB SNV/indel F1. Not in scope per the brief, but the design leaves a `--callers` param so it can be added as a second caller for comparison.
- *GATK gCNV for WGS* — needs a cohort/PoN; single-sample WGS is better served by CNVpytor + Manta.
- *Sentieon / DRAGEN* — commercial; out of scope.

---

## 2. Reference and resource bundle

Download once to `references/GRCh38/` (one-off `prepare_references.nf` workflow, ~2 h, dominated by the bwa-mem2 index build).

| Resource | Source | Notes |
|---|---|---|
| `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna` | NCBI | **No ALT contigs** (avoids MAPQ 0 in alt regions on a non-alt-aware pipeline), **with decoy** (hs38d1), chrM = rCRS. Same assembly used by DRAGEN and GIAB benchmarking. |
| `.fai`, `.dict`, bwa-mem2 index (`.bwt.2bit.64`, `.pac`, `.amb`, `.ann`, `.0123`) | generated | One index serves both aligners: bwa-mem3 reads the bwa-mem2 index as-is. Build it with `bwa-mem2 index` (or `bwa-mem3 index --emit-unpacked-ref`) so the `.0123` file exists for the bwa-mem2 fallback. |
| STR table (`ComposeSTRTableFile`) | generated | For DRAGEN mode. |
| dbSNP, Mills_and_1000G_gold_standard.indels, Homo_sapiens_assembly38.known_indels | GATK hg38 bundle (`gs://gcp-public-data--broad-references/hg38/v0`) | BQSR known sites. |
| `wgs_calling_regions.hg38.interval_list` | GATK bundle | WGS calling scatter base. |
| Exome capture BED(s) per kit | vendor (Twist, Agilent, IDT …) | Lift to GRCh38 if needed; padded ±100 bp for calling. |
| chrM resources (`chrM.fasta`, shifted-by-8000 fasta, `ShiftBack.chain`, blacklist, NuMT control regions) | GATK bundle `hg38/v0/chrM/` | Mito pipeline. |
| Haplogrep 3 / haplocheck tree data | bundled in containers | |
| GIAB HG002 WGS (30×) + HG002 WES FASTQ; GIAB v4.2.1 small-variant truth + bed; GIAB HG002 SV Tier 1 truth | NIST / GIAB FTP | Validation (§5). |

---

## 3. Pipeline architecture

```
samplesheet.csv
  sample,fastq_1,fastq_2,seq_type(wgs|wes),capture_bed,sex,lane   (multi-lane rows merged after alignment)
        │
        ▼
┌──────────────── QC_ALIGN (per lane) ─────────────────────────────┐
│ fastp ─▶ bwa-mem2 mem -Y -K 1e8 -t 64 │ samtools fixmate -m       │
│          │ samtools sort -@16 -m 2G    │ samtools markdup -d 2500 -S│  ──▶ CRAM + .crai
└──────────────────────────────────────────────────────────────────┘
        │  (lanes merged with samtools merge if >1)
        ├──────────────▶ ALIGN_QC: mosdepth, samtools stats, CollectWgs/HsMetrics,
        │                VerifyBamID2, somalier (sex) ──▶ QC gates
        ▼
┌──────────────── BQSR (scatter N shards, gather) ─────────────────┐
│ BaseRecalibrator ×N ─▶ GatherBQSRReports ─▶ ApplyBQSR ×N ─▶ merge│  ──▶ recal CRAM
└──────────────────────────────────────────────────────────────────┘
        │
        ├─────────────────────────────┬────────────────────────────┐
        ▼                             ▼                            ▼
GERMLINE_SNV                    MITO                          SV_CNV
 CalibrateDragstrModel           PrintReads -L chrM            wgs: Manta + CNVpytor
 HaplotypeCaller -ERC GVCF ×N    RevertSam → realign chrM      wes: Manta --exome
   --dragen-mode                 + shifted chrM                      + gCNV (case mode,
 MergeVcfs → GenotypeGVCFs       Mutect2 --mitochondria-mode ×2        requires PoN)
 VariantFiltration (DRAGEN)      LiftoverVcf → Merge → Filter   AnnotSV (optional)
 bcftools stats                  haplocheck, Haplogrep3
        │                             │                            │
        └─────────────────────────────┴────────────────────────────┘
                                      ▼
                                  MultiQC
```

**Outputs per sample** (`results/<sample>/`):
`alignment/*.cram(.crai)`, `qc/` (fastp, mosdepth, metrics, contamination, sex), `variants/<sample>.g.vcf.gz`, `variants/<sample>.filtered.vcf.gz`, `mito/<sample>.chrM.vcf.gz + haplogroup.txt + contamination.txt`, `sv_cnv/manta.*.vcf.gz`, `sv_cnv/cnvpytor.*.vcf`/`gcnv.*.vcf.gz`, `multiqc/`.

**Repository layout** (nf-core conventions, so `nf-core modules install` works):

```
main.nf
nextflow.config               # params, profiles (docker, apptainer, test, test_wes)
nextflow_schema.json          # param validation
conf/base.config              # resource labels (process_low/medium/high/high_memory)
conf/modules.config           # per-tool args + container versions (the single place to pin)
conf/test*.config             # chr20-subset CI configs
modules/nf-core/…             # installed nf-core modules (fastp, bwamem2, samtools, gatk4, manta, mosdepth, multiqc, …)
modules/local/…               # cnvpytor, haplocheck, haplogrep3, gatk mito helpers
subworkflows/local/           # qc_align, bqsr, germline_snv, mito, sv_cnv_wgs, sv_cnv_wes, pon_build
workflows/prepare_references.nf
assets/                       # samplesheet schema, multiqc_config.yml, interval lists
bin/                          # small helper scripts
tests/                        # nf-test
docs/
```

---

## 4. Speed design for one 80-core node

The goal is to keep all 80 cores busy for the whole run, not to make any one step fast in isolation.

1. **Stream alignment → sort → markdup**: one process, no SAM/BAM round-trip. bwa-mem2 `-t 64`, sort `-@ 16 -m 2G` (≈32 GB), markdup inline. This single step is ~40 % of WGS runtime.
2. **Scatter everything that GATK runs single-threaded**: BQSR and HaplotypeCaller over **20 shards** (`SplitIntervals --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION` on the WGS calling regions / padded capture BED). HC: 20 shards × 4 `--native-pair-hmm-threads` = 80 cores. Use `--smith-waterman FASTEST_AVAILABLE` and `--pair-hmm-implementation AVX_LOGLESS_CACHING_OMP`.
3. **Run Mito, Manta and CNVpytor concurrently** with HaplotypeCaller — Nextflow does this automatically; cap Manta at 16 threads so it doesn't starve HC.
4. **CRAM everywhere**, bgzip+tabix VCFs; Nextflow `cleanup = true` so the work dir doesn't eat the 1.1 TB (a 30× WGS run otherwise leaves ~300–400 GB of intermediates).
5. **Work dir on `/home`** (1.1 TB), *not* `/tmp` — `/tmp` is a 63 GB RAM-backed tmpfs and will fill on WGS.
6. **Multi-sample runs**: Nextflow's executor queue overlaps the alignment of sample *n+1* with the calling of sample *n*; `executor.cpus = 80`, `executor.memory = 115.GB`.

7. **Aligner: bwa-mem3 in `--compat=bwa-mem2` mode.** On multi-sample runs load the index once with `bwa-mem3 shm` (drop it with `shm -d` before any re-index — there is no staleness check). Alignment was the single largest step with bwa-mem2; bwa-mem3 should roughly halve the end-to-end WGS time if its benchmarks hold.

**Expected wall-clock (rough, 30× WGS, one sample, this node):** fastp+align+markdup ~1.5–2 h with bwa-mem2, **~30–45 min with bwa-mem3** (sort/markdup then dominate) · BQSR ~30 min · HC+genotype ~45–60 min · Mito ~10 min · Manta ~1 h (parallel) · CNVpytor ~30 min (parallel) ⇒ **≈ 2.5–3 h per WGS with bwa-mem3 (≈ 3.5–4.5 h with bwa-mem2)**, **≈ 20–45 min per WES**. Measured numbers replace these in Phase 6.

---

## 5. Accuracy and validation plan

- **Aligner head-to-head:** run HG002 WGS with `--aligner bwamem3` (compat mode) and `--aligner bwamem2`; alignment records must be byte-identical after stripping `MQ:i`/`HN:i`/`@PG` (`samtools view | cut` + `cmp`), and hap.py F1 must match to 4 decimals. Repeat on every bwa-mem3 version bump. Only then evaluate `--fast` and plain (non-compat) mode as a separate experiment.
- **Small variants:** run HG002 WGS (and HG002 WES) through the pipeline; compare to GIAB v4.2.1 with `hap.py` (or `rtg vcfeval`) inside the high-confidence BED, stratified (GIAB stratifications: homopolymers, segdups, MHC, …). Acceptance: SNV F1 ≥ 0.995, indel F1 ≥ 0.99 for WGS at 30× in the GIAB high-confidence region — typical for bwa-mem2 + HC DRAGEN mode.
- **SV/CNV:** `truvari bench` against GIAB HG002 SV Tier 1 (deletions/insertions ≥ 50 bp) for Manta; CNVpytor checked on known large HG002 CNVs and on an in-house sample with array/MLPA-confirmed CNVs when available.
- **Mito:** HG002 haplogroup must be reproduced by Haplogrep 3; haplocheck contamination < 1 %; heteroplasmy calls stable across two independent runs.
- **Regression:** a chr20 + chrM subset of HG002 (≈ 5 min end-to-end) runs in `nf-test` on every change (`-profile test`, `-profile test_wes`).
- **QC gates that fail a sample:** VerifyBamID2 FREEMIX > 0.03, sex mismatch (somalier vs. samplesheet), WGS mean coverage < 25× / WES ≥ 20× on < 90 % of target bases, duplication > 25 %.

---

## 6. Decisions to be aware of

| Decision | Choice | Trade-off |
|---|---|---|
| BQSR | **On by default**, `--skip_bqsr` available | Modern binned Illumina qualities make BQSR's benefit small, and the DRAGEN-GATK workflow omits it; it still costs ~30 min/WGS. Kept on as requested; measure its F1 effect in Phase 6 and decide. |
| Aligner | **bwa-mem3** in `--compat=bwa-mem2`, pinned; bwa-mem2 fallback | ~5× faster alignment and real bug fixes, at the cost of depending on a pre-1.0 project with breaking changes between minor releases and no published GIAB benchmark. Compat mode makes the output provably identical to bwa-mem2, so the residual risk is crashes/bugs, which the fallback covers. Never auto-bump the version. |
| Reference | GRCh38 **no-alt + decoy** | Loses HLA/alt-haplotype resolution, gains correct MAPQ in ~150 Mb of alt-shadowed regions without alt-aware post-processing. |
| Duplicate marking | samtools markdup | Fastest; results equivalent to Picard. Choose `MarkDuplicatesSpark` if a downstream consumer insists on Picard metrics format. |
| WES CNV | gCNV requires a PoN | No honest single-sample exome CNV without a reference set; the plan builds one from the lab's own data per kit. |
| Containers | Docker (needs group fix, see §7) with Apptainer as fallback profile | Both profiles ship so the pipeline runs on an HPC later. |

---

## 7. Implementation phases

| Phase | Deliverable | Est. |
|---|---|---|
| **0 – Environment** | Add user to `docker` group (`sudo usermod -aG docker labs`, re-login) *or* install Apptainer. Install Java 17 + Nextflow via micromamba (`micromamba create -n nf nextflow openjdk=17 nf-test`). Scaffold nf-core-style repo. | ½ day |
| **1 – References** | `prepare_references.nf`: download GRCh38 no-alt+decoy, build `.fai/.dict`, bwa-mem2 index, STR table, fetch GATK bundle, chrM resources, HG002 chr20 test subset. | ½ day + ~2 h compute |
| **2 – QC + alignment** | fastp → bwa-mem3 (`--compat=bwa-mem2`) / bwa-mem2 → samtools streaming → CRAM; `--aligner` switch; lane merge; ALIGN_QC subworkflow; MultiQC. Pass `-profile test` with both aligners. | 1–2 days |
| **3 – BQSR + germline SNV** | Scatter/gather BQSR; DRAGstr; HC GVCF (DRAGEN mode); GenotypeGVCFs; filtering; bcftools stats. WES interval handling. | 1–2 days |
| **4 – Mito** | Full GATK mito subworkflow + haplocheck + Haplogrep 3. | 1 day |
| **5 – SV/CNV** | Manta (wgs/exome), CNVpytor, gCNV case mode + `pon_build` workflow, optional AnnotSV. | 2 days |
| **6 – Validation + tuning** | HG002 WGS/WES full runs, hap.py / truvari reports, `-with-report` profiling, shard-count and thread tuning, finalize resource configs, measured runtimes in README. | 2–3 days |
| **7 – Hardening** | nf-test suite, `nextflow_schema.json`, README/usage docs, CHANGELOG, release `v1.0.0`. | 1 day |

Phases 2–5 are independent of each other after Phase 1 and can be developed in parallel branches.

---

## 8. Open questions for the lab

1. Which exome capture kit(s)? (Needed for the BED and the PoN.) How many WES samples already exist on that kit?
2. Sequencer / flow cell type (affects `markdup -d` optical distance: 2500 for patterned NovaSeq/NextSeq, 100 for HiSeq).
3. Is downstream annotation (VEP / snpEff, ClinVar, gnomAD) in scope for this pipeline, or handled elsewhere?
4. Retention policy: keep recalibrated CRAM only, or also pre-BQSR CRAM?
5. Will results feed a clinical report? If so, QC-gate thresholds in §5 need sign-off.
