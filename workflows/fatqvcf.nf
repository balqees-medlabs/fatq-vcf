/*
 * Main sample-level workflow. Subworkflows are added phase by phase (see docs/PIPELINE_PLAN.md §7).
 */

include { parseSamplesheet } from '../lib/samplesheet'
include { QC_ALIGN         } from '../subworkflows/local/qc_align'
include { MULTIQC          } from '../modules/local/multiqc/main'

workflow FATQVCF {
    take:
    samplesheet   // path: CSV

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    // ---- References ------------------------------------------------------------
    fasta          = file(params.fasta,          checkIfExists: true)
    fai            = file(params.fasta_fai,      checkIfExists: true)
    dict           = file(params.fasta_dict,     checkIfExists: true)
    bwa_index      = Channel.fromPath("${params.bwa_index}/*.{0123,amb,ann,bwt.2bit.64,pac}", checkIfExists: true).collect()
    somalier_sites = file(params.somalier_sites, checkIfExists: true)

    // ---- Samplesheet: one row per (sample, lane) ---------------------------------
    ch_reads = parseSamplesheet(samplesheet)
    ch_reads.count().subscribe { n -> log.info "Parsed ${n} FASTQ lane(s) from ${samplesheet}" }

    // ---- Phase 2: QC + alignment + alignment QC ----------------------------------
    QC_ALIGN( ch_reads, fasta, fai, dict, bwa_index, somalier_sites )
    ch_cram     = QC_ALIGN.out.cram
    ch_versions = ch_versions.mix(QC_ALIGN.out.versions)
    ch_multiqc  = ch_multiqc.mix(QC_ALIGN.out.multiqc)

    // Phase 3: BQSR + GERMLINE_SNV → Phase 4: MITO → Phase 5: SV_CNV

    // ---- Report -----------------------------------------------------------------
    ch_versions_yml = ch_versions.unique().collectFile(name: 'software_versions.yml', sort: true)
    MULTIQC( ch_multiqc.collect(), file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true), ch_versions_yml )

    emit:
    cram     = ch_cram
    multiqc  = MULTIQC.out.report
    versions = ch_versions
}
