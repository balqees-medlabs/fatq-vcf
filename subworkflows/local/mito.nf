/*
 * MITO: GATK mitochondria pipeline. chrM reads are re-aligned to both the standard and the
 * shifted-by-8000 chrM reference; Mutect2 --mitochondria-mode calls the non-control region on the
 * standard alignment and the control region on the shifted alignment; the shifted calls are lifted
 * back, merged, filtered, then haplocheck (contamination) + Haplogrep 3 (haplogroup) are run.
 */

include { MITO_SUBSET_REVERT } from '../../modules/local/mito/subset_revert'
include { MITO_ALIGN as MITO_ALIGN_PLAIN   } from '../../modules/local/mito/align'
include { MITO_ALIGN as MITO_ALIGN_SHIFTED } from '../../modules/local/mito/align'
include { MITO_MUTECT2 as MITO_MUTECT2_PLAIN   } from '../../modules/local/mito/mutect2'
include { MITO_MUTECT2 as MITO_MUTECT2_SHIFTED } from '../../modules/local/mito/mutect2'
include { MITO_MERGE_FILTER  } from '../../modules/local/mito/merge_filter'
include { HAPLOCHECK         } from '../../modules/local/mito/haplocheck'
include { HAPLOGREP3         } from '../../modules/local/mito/haplogrep3'

workflow MITO {
    take:
    ch_cram          // [ meta, cram, crai ]  (analysis-ready)
    fasta
    fai
    dict
    mt_plain         // [ fasta, fai, dict, [index files] ]   standard chrM
    mt_shifted       // [ fasta, fai, dict, [index files] ]   shifted-by-8000 chrM
    mt_non_control   // interval_list (standard coords)
    mt_control_shift // interval_list (shifted coords)
    shiftback_chain
    blacklist

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    // ---- chrM reads → unmapped BAM ---------------------------------------------
    MITO_SUBSET_REVERT( ch_cram, fasta, fai, dict )
    ch_ubam     = MITO_SUBSET_REVERT.out.ubam       // [ meta, ubam ]
    ch_versions = ch_versions.mix(MITO_SUBSET_REVERT.out.versions.first())

    // ---- Re-align to standard + shifted chrM -----------------------------------
    ch_plain_in = ch_ubam.map { meta, ubam ->
        [ meta, ubam, 'plain', mt_plain[0], mt_plain[1], mt_plain[2], mt_plain[3] ] }
    ch_shift_in = ch_ubam.map { meta, ubam ->
        [ meta, ubam, 'shifted', mt_shifted[0], mt_shifted[1], mt_shifted[2], mt_shifted[3] ] }

    MITO_ALIGN_PLAIN( ch_plain_in )
    MITO_ALIGN_SHIFTED( ch_shift_in )
    ch_versions = ch_versions.mix(MITO_ALIGN_PLAIN.out.versions.first())

    // ---- Mutect2 over the matching region --------------------------------------
    MITO_MUTECT2_PLAIN(
        MITO_ALIGN_PLAIN.out.bam.map { meta, tag, bam, bai ->
            [ meta, tag, bam, bai, mt_non_control, mt_plain[0], mt_plain[1], mt_plain[2] ] } )
    MITO_MUTECT2_SHIFTED(
        MITO_ALIGN_SHIFTED.out.bam.map { meta, tag, bam, bai ->
            [ meta, tag, bam, bai, mt_control_shift, mt_shifted[0], mt_shifted[1], mt_shifted[2] ] } )
    ch_versions = ch_versions.mix(MITO_MUTECT2_PLAIN.out.versions.first())

    // ---- Combine plain + shifted calls per sample ------------------------------
    ch_merge_in = MITO_MUTECT2_PLAIN.out.vcf
        .map { meta, tag, vcf, tbi, stats -> [ meta, vcf, tbi, stats ] }
        .join( MITO_MUTECT2_SHIFTED.out.vcf.map { meta, tag, vcf, tbi, stats -> [ meta, vcf, tbi, stats ] } )

    MITO_MERGE_FILTER( ch_merge_in, mt_plain[0], mt_plain[1], mt_plain[2], shiftback_chain, blacklist )
    ch_mito_vcf = MITO_MERGE_FILTER.out.vcf          // [ meta, vcf, tbi ]
    ch_versions = ch_versions.mix(MITO_MERGE_FILTER.out.versions.first())

    // ---- Contamination + haplogroup --------------------------------------------
    HAPLOCHECK( ch_mito_vcf )
    HAPLOGREP3( ch_mito_vcf )
    ch_versions = ch_versions.mix(HAPLOCHECK.out.versions.first(), HAPLOGREP3.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(HAPLOCHECK.out.report.map { m, f -> f })

    emit:
    vcf        = ch_mito_vcf                          // [ meta, vcf, tbi ]
    haplocheck = HAPLOCHECK.out.report
    haplogroup = HAPLOGREP3.out.haplogroup
    multiqc    = ch_multiqc
    versions   = ch_versions
}
