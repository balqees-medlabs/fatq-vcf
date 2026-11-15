/*
 * PON_BUILD (Phase 7 helper): build a GATK germline-CNV panel of normals for a WES capture kit.
 *
 *   nextflow run main.nf -profile docker --pon_build \
 *     --pon_input pon_samplesheet.csv --pon_bed kit.hg38.bed --pon_name twist2
 *
 * pon_samplesheet.csv: one column `cram` listing ≥ 20–30 analysis-ready CRAMs from the same kit.
 * Produces: preprocessed intervals, the ploidy model and the gCNV model — the inputs that
 * --gcnv_ploidy_model_dir / --gcnv_model_dir / --gcnv_intervals expect for case-mode calling.
 */

include { GATK4_PREPROCESSINTERVALS      } from '../modules/local/gatk4/preprocessintervals'
include { GATK4_COLLECTREADCOUNTS as PON_COUNTS } from '../modules/local/gatk4/collectreadcounts'
include { GATK4_FILTERINTERVALS         } from '../modules/local/gatk4/filterintervals'
include { GATK4_CONTIGPLOIDY_COHORT     } from '../modules/local/gatk4/contigploidy_cohort'
include { GATK4_GERMLINECNV_COHORT      } from '../modules/local/gatk4/gcnv_cohort'

workflow PON_BUILD {
    main:
    fasta = file(params.fasta,      checkIfExists: true)
    fai   = file(params.fasta_fai,  checkIfExists: true)
    dict  = file(params.fasta_dict, checkIfExists: true)
    bed   = file(params.pon_bed,    checkIfExists: true)

    if (!params.pon_input) error "ERROR: --pon_input samplesheet (column 'cram') is required with --pon_build"

    // Analysis-ready CRAMs from the cohort
    ch_crams = Channel.fromPath(params.pon_input)
        .splitCsv(header: true)
        .map { row -> def c = file(row.cram, checkIfExists: true)
                      [ [ id: c.simpleName ], c, file("${c}.crai", checkIfExists: true) ] }

    GATK4_PREPROCESSINTERVALS( bed, fasta, fai, dict )
    PON_COUNTS( ch_crams, fasta, fai, dict, GATK4_PREPROCESSINTERVALS.out.intervals )

    ch_counts = PON_COUNTS.out.counts.map { meta, h -> h }.collect()
    GATK4_FILTERINTERVALS( GATK4_PREPROCESSINTERVALS.out.intervals, ch_counts, fasta, fai, dict )
    GATK4_CONTIGPLOIDY_COHORT( ch_counts, GATK4_FILTERINTERVALS.out.intervals )
    GATK4_GERMLINECNV_COHORT( ch_counts, GATK4_FILTERINTERVALS.out.intervals, GATK4_CONTIGPLOIDY_COHORT.out.calls )
}
