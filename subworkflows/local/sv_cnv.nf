/*
 * SV_CNV: structural + copy-number calling.
 *   WGS → Manta (SV) + CNVpytor (read-depth CNV)
 *   WES → Manta --exome (SV) + GATK gCNV case mode (CNV; only when a PoN model is provided)
 * Optional AnnotSV annotation of the Manta diploid SVs when the annotation bundle is given.
 */

include { MANTA_GERMLINE          } from '../../modules/local/manta/germline'
include { CNVPYTOR                } from '../../modules/local/cnvpytor/main'
include { GATK4_COLLECTREADCOUNTS } from '../../modules/local/gatk4/collectreadcounts'
include { GATK4_GERMLINECNV_CASE  } from '../../modules/local/gatk4/gcnv_case'
include { ANNOTSV                 } from '../../modules/local/annotsv/main'

workflow SV_CNV {
    take:
    ch_cram          // [ meta, cram, crai ]
    ch_snv_vcf       // [ meta, vcf, tbi ]  (for CNVpytor BAF)
    fasta
    fai
    dict

    main:
    ch_versions = Channel.empty()

    // ---- Manta (all samples; --exome added for WES via conf/modules.config) -----
    MANTA_GERMLINE( ch_cram.map { meta, cram, crai -> [ meta, cram, crai, [], [] ] }, fasta, fai )
    ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions.first())

    // ---- WGS: CNVpytor ----------------------------------------------------------
    ch_wgs = ch_cram.filter { meta, cram, crai -> meta.seq_type == 'wgs' }
        .join( ch_snv_vcf, remainder: true )
        .map { meta, cram, crai, vcf, tbi -> [ meta, cram, crai, vcf ?: [], tbi ?: [] ] }
    CNVPYTOR( ch_wgs, fasta, fai )
    ch_versions = ch_versions.mix(CNVPYTOR.out.versions.first())

    // ---- WES: GATK gCNV case mode (only if a PoN model is configured) -----------
    ch_gcnv_vcf = Channel.empty()
    if (params.gcnv_model_dir && params.gcnv_ploidy_model_dir) {
        ch_wes = ch_cram.filter { meta, cram, crai -> meta.seq_type == 'wes' }
        GATK4_COLLECTREADCOUNTS( ch_wes, fasta, fai, dict, file(params.gcnv_intervals, checkIfExists: true) )
        GATK4_GERMLINECNV_CASE(
            GATK4_COLLECTREADCOUNTS.out.counts,
            file(params.gcnv_ploidy_model_dir, checkIfExists: true),
            file(params.gcnv_model_dir, checkIfExists: true),
            dict )
        ch_gcnv_vcf = GATK4_GERMLINECNV_CASE.out.vcf
        ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions.first(), GATK4_GERMLINECNV_CASE.out.versions.first())
    }

    // ---- Optional AnnotSV on the Manta diploid SVs ------------------------------
    if (params.annotsv_dir) {
        ANNOTSV( MANTA_GERMLINE.out.diploid_sv.map { meta, vcf, tbi -> [ meta, vcf ] },
                 file(params.annotsv_dir, checkIfExists: true) )
        ch_versions = ch_versions.mix(ANNOTSV.out.versions.first())
    }

    emit:
    manta_sv  = MANTA_GERMLINE.out.diploid_sv
    cnvpytor  = CNVPYTOR.out.vcf
    gcnv      = ch_gcnv_vcf
    versions  = ch_versions
}
