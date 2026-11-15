/*
 * VALIDATE (Phase 6): benchmark pipeline outputs against GIAB truth.
 *
 *   nextflow run main.nf -profile docker --validate \
 *     --validate_vcf results/HG002/variants/HG002.filtered.vcf.gz \
 *     --truth_vcf  references/GRCh38/validation/giab_v4.2.1/HG002_...benchmark.vcf.gz \
 *     --truth_bed  references/GRCh38/validation/giab_v4.2.1/HG002_...benchmark_noinconsistent.bed
 *
 * Optional SV benchmark (--sv_vcf / --sv_truth_vcf / --sv_truth_bed) runs truvari.
 */

include { HAPPY         } from '../modules/local/validation/happy'
include { TRUVARI_BENCH } from '../modules/local/validation/truvari'

workflow VALIDATE {
    main:
    fasta = file(params.fasta,     checkIfExists: true)
    fai   = file(params.fasta_fai, checkIfExists: true)

    if (!params.validate_vcf) error "ERROR: --validate_vcf is required with --validate"

    // ---- Small-variant benchmark (hap.py) --------------------------------------
    def vcf = file(params.validate_vcf, checkIfExists: true)
    def tbi = file("${params.validate_vcf}.tbi", checkIfExists: true)
    def meta = [ id: params.validate_id ?: vcf.simpleName ]

    ch_query = Channel.of([ meta, vcf, tbi ])
    ch_truth = Channel.of([
        file(params.truth_vcf, checkIfExists: true),
        file("${params.truth_vcf}.tbi", checkIfExists: true),
        file(params.truth_bed, checkIfExists: true) ]).collect()

    HAPPY( ch_query, ch_truth, fasta, fai )
    HAPPY.out.summary.subscribe { m, csv -> log.info "hap.py summary for ${m.id}: ${csv}" }

    // ---- Optional SV benchmark (truvari) ---------------------------------------
    if (params.sv_vcf && params.sv_truth_vcf && params.sv_truth_bed) {
        def sv  = file(params.sv_vcf, checkIfExists: true)
        def svt = file("${params.sv_vcf}.tbi", checkIfExists: true)
        ch_sv_query = Channel.of([ meta, sv, svt ])
        ch_sv_truth = Channel.of([
            file(params.sv_truth_vcf, checkIfExists: true),
            file("${params.sv_truth_vcf}.tbi", checkIfExists: true),
            file(params.sv_truth_bed, checkIfExists: true) ]).collect()
        TRUVARI_BENCH( ch_sv_query, ch_sv_truth, fasta, fai )
    }
}
