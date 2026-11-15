/*
 * QC_ALIGN: per-lane fastp→aligner→sort, per-sample lane merge + markdup → CRAM,
 *           then alignment QC (samtools stats, mosdepth, VerifyBamID2, somalier, HsMetrics for WES)
 *           and sample-level QC gates.
 */

include { FASTP_BWAMEM           } from '../../modules/local/align/fastp_bwamem'
include { SAMTOOLS_MARKDUP       } from '../../modules/local/samtools/markdup'
include { SAMTOOLS_STATS         } from '../../modules/local/samtools/stats'
include { MOSDEPTH               } from '../../modules/local/mosdepth/main'
include { VERIFYBAMID2           } from '../../modules/local/verifybamid2/main'
include { SOMALIER_EXTRACT       } from '../../modules/local/somalier/main'
include { SOMALIER_RELATE        } from '../../modules/local/somalier/main'
include { PICARD_COLLECTHSMETRICS } from '../../modules/local/picard/collecthsmetrics'
include { QC_GATES               } from '../../modules/local/qc_gates/main'

workflow QC_ALIGN {
    take:
    ch_reads        // [ meta(lane-level), [fq1, fq2] ]
    fasta
    fai
    dict
    bwa_index       // list of index files
    somalier_sites

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    // ---- Align each lane -----------------------------------------------------
    FASTP_BWAMEM( ch_reads, fasta, fai, bwa_index )
    ch_versions = ch_versions.mix(FASTP_BWAMEM.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(FASTP_BWAMEM.out.fastp_json.map { meta, f -> f })

    // ---- Group lanes per sample; drop lane from meta ---------------------------
    ch_sample_bams = FASTP_BWAMEM.out.bam
        .map { meta, bam -> [ meta.subMap(['id', 'sample', 'seq_type', 'capture_bed', 'sex', 'single_end']), bam ] }
        .groupTuple()

    SAMTOOLS_MARKDUP( ch_sample_bams, fasta, fai )
    ch_cram     = SAMTOOLS_MARKDUP.out.cram                       // [ meta, cram, crai ]
    ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(SAMTOOLS_MARKDUP.out.stats.map { meta, f -> f })

    // ---- Alignment QC ----------------------------------------------------------
    SAMTOOLS_STATS( ch_cram, fasta, fai )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(
        SAMTOOLS_STATS.out.stats.map { m, f -> f },
        SAMTOOLS_STATS.out.flagstat.map { m, f -> f },
        SAMTOOLS_STATS.out.idxstats.map { m, f -> f })

    // mosdepth: WES by capture BED, WGS by 500 bp windows (bed slot empty)
    MOSDEPTH( ch_cram.map { meta, cram, crai -> [ meta, cram, crai, meta.seq_type == 'wes' ? meta.capture_bed : [] ] }, fasta, fai )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(
        MOSDEPTH.out.summary.map { m, f -> f },
        MOSDEPTH.out.global_dist.map { m, f -> f },
        MOSDEPTH.out.region_dist.map { m, f -> f })

    VERIFYBAMID2( ch_cram, fasta, fai, params.verifybamid_svd )
    ch_versions = ch_versions.mix(VERIFYBAMID2.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(VERIFYBAMID2.out.selfsm.map { m, f -> f })

    SOMALIER_EXTRACT( ch_cram, fasta, fai, somalier_sites )
    SOMALIER_RELATE( SOMALIER_EXTRACT.out.extract.map { m, f -> f }.collect() )
    ch_versions = ch_versions.mix(SOMALIER_EXTRACT.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(SOMALIER_RELATE.out.samples, SOMALIER_RELATE.out.pairs)

    // WES only: Picard HsMetrics
    PICARD_COLLECTHSMETRICS(
        ch_cram.filter { meta, cram, crai -> meta.seq_type == 'wes' }
               .map { meta, cram, crai -> [ meta, cram, crai, meta.capture_bed ] },
        fasta, fai, dict )
    ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(PICARD_COLLECTHSMETRICS.out.metrics.map { m, f -> f })

    // ---- QC gates ----------------------------------------------------------------
    ch_gate_in = MOSDEPTH.out.summary
        .join( MOSDEPTH.out.thresholds.map { m, f -> [m, f] }, remainder: true )
        .map { meta, summ, thr -> [ meta, summ, thr ?: [] ] }
        .join( SAMTOOLS_MARKDUP.out.stats )
        .join( VERIFYBAMID2.out.selfsm )

    QC_GATES( ch_gate_in, SOMALIER_RELATE.out.samples.first() )
    ch_versions = ch_versions.mix(QC_GATES.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(QC_GATES.out.mqc)

    emit:
    cram       = ch_cram                    // [ meta, cram, crai ]
    qc_summary = QC_GATES.out.summary       // [ meta, tsv ]
    multiqc    = ch_multiqc
    versions   = ch_versions
}
