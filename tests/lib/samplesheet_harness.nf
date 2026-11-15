// Test harness: expose parseSamplesheet's channel as a named workflow output for nf-test.
include { parseSamplesheet } from '../../lib/samplesheet'

workflow SAMPLESHEET_HARNESS {
    take:
    sheet

    main:
    rows = parseSamplesheet(sheet)

    emit:
    rows
}
