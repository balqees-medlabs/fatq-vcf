/*
 * BQSR: BaseRecalibrator (scatter over pre-split shards) → GatherBQSRReports →
 *       ApplyBQSR (scatter) → merge recalibrated CRAM.
 *
 * Shards come from PREPARE_INTERVALS (shared with GERMLINE_SNV).
 */

include { GATK4_BASERECALIBRATOR  } from '../../modules/local/gatk4/baserecalibrator'
include { GATK4_GATHERBQSRREPORTS } from '../../modules/local/gatk4/gatherbqsrreports'
include { GATK4_APPLYBQSR         } from '../../modules/local/gatk4/applybqsr'
include { SAMTOOLS_MERGE_CRAM     } from '../../modules/local/samtools/merge_cram'

workflow BQSR {
    take:
    ch_cram          // [ meta, cram, crai ]
    ch_shards        // [ meta, shard ]
    fasta
    fai
    dict
    dbsnp            // [ vcf, tbi ]
    mills            // [ vcf, tbi ]
    known_indels     // [ vcf, tbi ]

    main:
    ch_versions = Channel.empty()

    // ---- BaseRecalibrator per shard → gather -----------------------------------
    ch_br_in = ch_cram.combine(ch_shards, by: 0)      // [ meta, cram, crai, shard ]
    GATK4_BASERECALIBRATOR( ch_br_in, fasta, fai, dict, dbsnp, mills, known_indels )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first())

    GATK4_GATHERBQSRREPORTS( GATK4_BASERECALIBRATOR.out.table.groupTuple() )
    ch_table = GATK4_GATHERBQSRREPORTS.out.table       // [ meta, table ]

    // ---- ApplyBQSR per shard → merge -------------------------------------------
    ch_ab_in = ch_cram.combine(ch_table, by: 0)        // [ meta, cram, crai, table ]
        .combine(ch_shards, by: 0)                     // [ meta, cram, crai, table, shard ]
    GATK4_APPLYBQSR( ch_ab_in, fasta, fai, dict )
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions.first())

    SAMTOOLS_MERGE_CRAM( GATK4_APPLYBQSR.out.cram.groupTuple(), fasta, fai )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_CRAM.out.versions.first())

    emit:
    cram     = SAMTOOLS_MERGE_CRAM.out.cram            // [ meta, cram, crai ]
    versions = ch_versions
}
