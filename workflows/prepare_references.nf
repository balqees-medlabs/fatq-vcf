/*
 * One-off reference bundle builder (docs/PIPELINE_PLAN.md §2).
 *
 *   nextflow run main.nf -profile docker --prepare_references [--test_data false] [--validation_data]
 *
 * Writes into params.refdir:
 *   GRCh38 no-alt + hs38d1 decoy FASTA, .fai, .dict, bwa-mem2/3 index, STR table
 *   gatk_bundle/   dbSNP (bgzipped), Mills, known_indels, WGS calling regions
 *   chrM/          Broad mito resources (+ indexes for chrM and shifted chrM)
 *   somalier/      sites.hg38.vcf.gz
 * and, with --test_data, a HG002 chr20+chrM ~30x FASTQ pair + GIAB chr20 truth into tests/data/.
 */

include { FETCH_GENOME                   } from '../modules/local/references/fetch_genome'
include { FETCH_URL                      } from '../modules/local/references/fetch_url'
include { SAMTOOLS_FAIDX                 } from '../modules/local/samtools/faidx'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/local/gatk4/createsequencedictionary'
include { BWAMEM2_INDEX                  } from '../modules/local/bwamem2/index'
include { GATK4_COMPOSESTRTABLEFILE      } from '../modules/local/gatk4/composestrtablefile'
include { BGZIP_TABIX                    } from '../modules/local/references/bgzip_tabix'
include { CHRM_INDEX                     } from '../modules/local/references/chrm_index'
include { GIAB_TEST_SUBSET               } from '../modules/local/references/giab_test_subset'
include { GIAB_TRUTH_SUBSET              } from '../modules/local/references/giab_truth_subset'
// Alias so the validation downloads get their own publishDir without colliding with FETCH_URL.
include { FETCH_URL as FETCH_URL_VALIDATION } from '../modules/local/references/fetch_url'
include { FETCH_URL as FETCH_URL_BED        } from '../modules/local/references/fetch_url'

workflow PREPARE_REFERENCES {

    main:
    def BROAD = 'https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0'
    def GIAB  = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp'

    // ---- Genome --------------------------------------------------------------
    FETCH_GENOME(
        'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz'
    )
    SAMTOOLS_FAIDX( FETCH_GENOME.out.fasta )
    GATK4_CREATESEQUENCEDICTIONARY( FETCH_GENOME.out.fasta )
    BWAMEM2_INDEX( FETCH_GENOME.out.fasta )
    GATK4_COMPOSESTRTABLEFILE( FETCH_GENOME.out.fasta, SAMTOOLS_FAIDX.out.fai, GATK4_CREATESEQUENCEDICTIONARY.out.dict )

    // ---- Flat downloads: [meta(subdir, name), url] ---------------------------
    ch_downloads = Channel.of(
        [ 'gatk_bundle', "${BROAD}/Homo_sapiens_assembly38.dbsnp138.vcf" ],
        [ 'gatk_bundle', "${BROAD}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" ],
        [ 'gatk_bundle', "${BROAD}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" ],
        [ 'gatk_bundle', "${BROAD}/Homo_sapiens_assembly38.known_indels.vcf.gz" ],
        [ 'gatk_bundle', "${BROAD}/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi" ],
        [ 'gatk_bundle', "${BROAD}/wgs_calling_regions.hg38.interval_list" ],
        [ 'chrM',        "${BROAD}/chrM/Homo_sapiens_assembly38.chrM.fasta" ],
        [ 'chrM',        "${BROAD}/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta" ],
        [ 'chrM',        "${BROAD}/chrM/ShiftBack.chain" ],
        [ 'chrM',        "${BROAD}/chrM/blacklist_sites.hg38.chrM.bed" ],
        [ 'chrM',        "${BROAD}/chrM/control_region_shifted.chrM.interval_list" ],
        [ 'chrM',        "${BROAD}/chrM/non_control_region.chrM.interval_list" ],
        [ 'somalier',    'https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz' ],
    ).map { subdir, url -> [ [ subdir: subdir, name: url.tokenize('/')[-1] ], url ] }

    FETCH_URL( ch_downloads )

    // dbSNP ships as an 11 GB plain VCF: bgzip + tabix it (GATK reads .vcf.gz natively).
    FETCH_URL.out.file
        .filter { meta, f -> meta.name == 'Homo_sapiens_assembly38.dbsnp138.vcf' }
        | BGZIP_TABIX

    // chrM + shifted chrM need .fai/.dict/bwa index for the mito realignment step.
    FETCH_URL.out.file
        .filter { meta, f -> meta.subdir == 'chrM' && meta.name.endsWith('.fasta') }
        | CHRM_INDEX

    // ---- Optional: HG002 chr20+chrM test subset + GIAB truth -----------------
    if (params.test_data) {
        GIAB_TEST_SUBSET(
            "${GIAB}/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam",
            FETCH_GENOME.out.fasta,
            SAMTOOLS_FAIDX.out.fai,
            GATK4_CREATESEQUENCEDICTIONARY.out.dict
        )
        // Full benchmark BED fetched on the host (bcftools container has no curl), then filtered.
        ch_truth_bed = Channel.of(
            [ [ subdir: 'validation/giab_v4.2.1', name: 'HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed' ],
              "${GIAB}/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" ] )
        FETCH_URL_BED( ch_truth_bed )
        GIAB_TRUTH_SUBSET(
            "${GIAB}/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
            FETCH_URL_BED.out.file.map { meta, f -> f }
        )
    }

    // ---- Optional: full HG002 30x WGS FASTQ for Phase 6 validation (52 GB) --
    if (params.validation_data) {
        def BRAIN = 'https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x'
        Channel.of(
            [ 'validation/HG002_WGS', "${BRAIN}/HG002.novaseq.pcr-free.30x.R1.fastq.gz" ],
            [ 'validation/HG002_WGS', "${BRAIN}/HG002.novaseq.pcr-free.30x.R2.fastq.gz" ],
            [ 'validation/giab_v4.2.1', "${GIAB}/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" ],
            [ 'validation/giab_v4.2.1', "${GIAB}/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" ],
            [ 'validation/giab_v4.2.1', "${GIAB}/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" ],
        ).map { subdir, url -> [ [ subdir: subdir, name: url.tokenize('/')[-1] ], url ] }
        | FETCH_URL_VALIDATION
    }
}
