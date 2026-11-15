/*
 * Main sample-level workflow. Subworkflows are added phase by phase (see docs/PIPELINE_PLAN.md §7).
 */

include { parseSamplesheet   } from '../lib/samplesheet'
include { QC_ALIGN           } from '../subworkflows/local/qc_align'
include { PREPARE_INTERVALS  } from '../subworkflows/local/prepare_intervals'
include { BQSR               } from '../subworkflows/local/bqsr'
include { GERMLINE_SNV       } from '../subworkflows/local/germline_snv'
include { MITO               } from '../subworkflows/local/mito'
include { SV_CNV             } from '../subworkflows/local/sv_cnv'
include { MULTIQC            } from '../modules/local/multiqc/main'

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
    str_table      = file(params.str_table,      checkIfExists: true)

    wgs_calling_regions = file(params.wgs_calling_regions, checkIfExists: true)
    override_intervals  = params.intervals ? file(params.intervals, checkIfExists: true) : []

    dbsnp        = [ file(params.dbsnp,        checkIfExists: true), file(params.dbsnp_tbi,        checkIfExists: true) ]
    mills        = [ file(params.mills,        checkIfExists: true), file(params.mills_tbi,        checkIfExists: true) ]
    known_indels = [ file(params.known_indels, checkIfExists: true), file(params.known_indels_tbi, checkIfExists: true) ]

    // ---- Samplesheet: one row per (sample, lane) ---------------------------------
    ch_reads = parseSamplesheet(samplesheet)
    ch_reads.count().subscribe { n -> log.info "Parsed ${n} FASTQ lane(s) from ${samplesheet}" }

    // ---- Phase 2: QC + alignment + alignment QC ----------------------------------
    QC_ALIGN( ch_reads, fasta, fai, dict, bwa_index, somalier_sites )
    ch_markdup_cram = QC_ALIGN.out.cram
    ch_versions     = ch_versions.mix(QC_ALIGN.out.versions)
    ch_multiqc      = ch_multiqc.mix(QC_ALIGN.out.multiqc)

    // ---- Phase 3: BQSR + germline SNV/indel --------------------------------------
    // Split calling intervals once; the shards feed both BQSR and HaplotypeCaller.
    PREPARE_INTERVALS( ch_markdup_cram, fasta, fai, dict, wgs_calling_regions, override_intervals )
    ch_shards   = PREPARE_INTERVALS.out.shards
    ch_versions = ch_versions.mix(PREPARE_INTERVALS.out.versions)

    if (params.skip_bqsr) {
        ch_analysis_cram = ch_markdup_cram
    } else {
        BQSR( ch_markdup_cram, ch_shards, fasta, fai, dict, dbsnp, mills, known_indels )
        ch_analysis_cram = BQSR.out.cram
        ch_versions      = ch_versions.mix(BQSR.out.versions)
    }

    GERMLINE_SNV( ch_analysis_cram, ch_shards, fasta, fai, dict, str_table, dbsnp )
    ch_versions = ch_versions.mix(GERMLINE_SNV.out.versions)
    ch_multiqc  = ch_multiqc.mix(GERMLINE_SNV.out.multiqc)

    // ---- Phase 4: mitochondrial calling ------------------------------------------
    if (!params.skip_mito) {
        def mt_plain = [
            file(params.mt_fasta, checkIfExists: true),
            file("${params.mt_fasta}.fai", checkIfExists: true),
            file("${params.mt_dir}/Homo_sapiens_assembly38.chrM.dict", checkIfExists: true),
            file("${params.mt_fasta}.{0123,amb,ann,bwt.2bit.64,pac}", checkIfExists: true) ]
        def mt_shifted = [
            file(params.mt_shifted_fasta, checkIfExists: true),
            file("${params.mt_shifted_fasta}.fai", checkIfExists: true),
            file("${params.mt_dir}/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict", checkIfExists: true),
            file("${params.mt_shifted_fasta}.{0123,amb,ann,bwt.2bit.64,pac}", checkIfExists: true) ]
        MITO(
            ch_analysis_cram, fasta, fai, dict,
            mt_plain, mt_shifted,
            file(params.mt_non_control_region,     checkIfExists: true),
            file(params.mt_control_region_shifted, checkIfExists: true),
            file(params.mt_shiftback_chain,        checkIfExists: true),
            file(params.mt_blacklist,              checkIfExists: true) )
        ch_versions = ch_versions.mix(MITO.out.versions)
        ch_multiqc  = ch_multiqc.mix(MITO.out.multiqc)
    }

    // ---- Phase 5: SV / CNV -------------------------------------------------------
    if (!params.skip_sv_cnv) {
        SV_CNV( ch_analysis_cram, GERMLINE_SNV.out.vcf, fasta, fai, dict )
        ch_versions = ch_versions.mix(SV_CNV.out.versions)
    }

    // ---- Report -----------------------------------------------------------------
    ch_versions_yml = ch_versions.unique().collectFile(name: 'software_versions.yml', sort: true)
    MULTIQC( ch_multiqc.collect(), file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true), ch_versions_yml )

    emit:
    cram     = ch_analysis_cram
    vcf      = GERMLINE_SNV.out.vcf
    gvcf     = GERMLINE_SNV.out.gvcf
    multiqc  = MULTIQC.out.report
    versions = ch_versions
}
