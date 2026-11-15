/*
 * GERMLINE_SNV: CalibrateDragstrModel → HaplotypeCaller (GVCF + DRAGEN mode, scatter) →
 *               MergeVcfs (GVCF) → GenotypeGVCFs → VariantFiltration (DRAGEN hard filter) →
 *               bcftools stats.
 *
 * Reuses the shard channel produced by BQSR so intervals are split only once.
 */

include { GATK4_CALIBRATEDRAGSTRMODEL } from '../../modules/local/gatk4/calibratedragstrmodel'
include { GATK4_HAPLOTYPECALLER       } from '../../modules/local/gatk4/haplotypecaller'
include { GATK4_MERGEVCFS             } from '../../modules/local/gatk4/mergevcfs'
include { GATK4_GENOTYPEGVCFS         } from '../../modules/local/gatk4/genotypegvcfs'
include { GATK4_VARIANTFILTRATION     } from '../../modules/local/gatk4/variantfiltration'
include { BCFTOOLS_STATS              } from '../../modules/local/bcftools/stats'

workflow GERMLINE_SNV {
    take:
    ch_cram          // [ meta, cram, crai ]  (recalibrated)
    ch_shards        // [ meta, shard ]
    fasta
    fai
    dict
    str_table
    dbsnp            // [ vcf, tbi ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    // ---- DRAGstr model per sample ----------------------------------------------
    GATK4_CALIBRATEDRAGSTRMODEL( ch_cram, fasta, fai, dict, str_table )
    ch_versions = ch_versions.mix(GATK4_CALIBRATEDRAGSTRMODEL.out.versions.first())

    // ---- HaplotypeCaller per shard ---------------------------------------------
    ch_hc_in = ch_cram
        .combine(GATK4_CALIBRATEDRAGSTRMODEL.out.model, by: 0)   // [ meta, cram, crai, model ]
        .combine(ch_shards, by: 0)                               // [ meta, cram, crai, model, shard ]
    GATK4_HAPLOTYPECALLER( ch_hc_in, fasta, fai, dict, dbsnp )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    // ---- Gather GVCF shards (join vcf+tbi lists per sample) --------------------
    ch_gvcf_parts = GATK4_HAPLOTYPECALLER.out.gvcf
        .map { meta, vcf, tbi -> [ meta, vcf, tbi ] }
        .groupTuple()
    GATK4_MERGEVCFS( ch_gvcf_parts, dict )
    ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions.first())

    // ---- Genotype → filter → stats ---------------------------------------------
    GATK4_GENOTYPEGVCFS( GATK4_MERGEVCFS.out.vcf, fasta, fai, dict, dbsnp )
    GATK4_VARIANTFILTRATION( GATK4_GENOTYPEGVCFS.out.vcf, fasta, fai, dict )
    BCFTOOLS_STATS( GATK4_VARIANTFILTRATION.out.vcf )
    ch_versions = ch_versions.mix(
        GATK4_GENOTYPEGVCFS.out.versions.first(),
        GATK4_VARIANTFILTRATION.out.versions.first(),
        BCFTOOLS_STATS.out.versions.first() )
    ch_multiqc = ch_multiqc.mix(BCFTOOLS_STATS.out.stats.map { m, f -> f })

    emit:
    gvcf       = GATK4_MERGEVCFS.out.vcf                 // [ meta, g.vcf.gz, tbi ]
    vcf        = GATK4_VARIANTFILTRATION.out.vcf         // [ meta, filtered.vcf.gz, tbi ]
    multiqc    = ch_multiqc
    versions   = ch_versions
}
