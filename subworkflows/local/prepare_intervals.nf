/*
 * PREPARE_INTERVALS: per-sample calling intervals (WGS regions or padded WES capture BED) split
 * into N balanced shards. Run once; the shard channel feeds both BQSR and GERMLINE_SNV so the
 * genome is split a single time per sample.
 */

include { INTERVALS_PREPARE    } from '../../modules/local/gatk4/intervals_prepare'
include { GATK4_SPLITINTERVALS } from '../../modules/local/gatk4/splitintervals'

workflow PREPARE_INTERVALS {
    take:
    ch_cram          // [ meta, cram, crai ]
    fasta
    fai
    dict
    wgs_calling_regions
    override_intervals

    main:
    ch_versions = Channel.empty()

    INTERVALS_PREPARE(
        ch_cram.map { meta, cram, crai -> [ meta, meta.seq_type == 'wes' ? meta.capture_bed : [] ] },
        dict, wgs_calling_regions, override_intervals )
    GATK4_SPLITINTERVALS( INTERVALS_PREPARE.out.intervals, fasta, fai, dict )
    ch_versions = ch_versions.mix(INTERVALS_PREPARE.out.versions.first(), GATK4_SPLITINTERVALS.out.versions.first())

    ch_shards = GATK4_SPLITINTERVALS.out.shards
        .flatMap { meta, shards -> (shards instanceof List ? shards : [shards]).collect { [ meta, it ] } }

    emit:
    shards   = ch_shards     // [ meta, shard_interval_list ]
    versions = ch_versions
}
