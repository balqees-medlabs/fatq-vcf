/*
 * Samplesheet parsing + validation (plain Groovy; nf-schema validation arrives in Phase 7).
 *
 * Columns: sample,lane,fastq_1,fastq_2,seq_type,capture_bed,sex
 */

def parseSamplesheet(samplesheet) {
    def required = ['sample', 'fastq_1', 'seq_type']
    def seen     = [] as Set

    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            required.each { col ->
                if (!row.containsKey(col) || !row[col]) error "Samplesheet: column '${col}' missing or empty in row ${row}"
            }
            def sample   = row.sample
            def lane     = row.lane ?: 'L001'
            def seq_type = row.seq_type.toLowerCase()
            if (!(seq_type in ['wgs', 'wes'])) error "Samplesheet: seq_type must be wgs or wes (sample ${sample}: '${row.seq_type}')"

            def sex = (row.sex ?: 'unknown').toLowerCase()
            if (!(sex in ['male', 'female', 'unknown'])) error "Samplesheet: sex must be male|female|unknown (sample ${sample}: '${row.sex}')"

            def key = "${sample}\t${lane}"
            if (key in seen) error "Samplesheet: duplicate sample/lane pair ${sample}/${lane}"
            seen << key

            def fq1 = file(row.fastq_1, checkIfExists: true)
            def fq2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
            def bed = row.capture_bed ? file(row.capture_bed, checkIfExists: true) : null
            if (seq_type == 'wes' && !bed) error "Samplesheet: wes sample ${sample} requires capture_bed"

            def meta = [
                id         : sample,
                sample     : sample,
                lane       : lane,
                seq_type   : seq_type,
                capture_bed: bed,
                sex        : sex,
                single_end : fq2 == null,
            ]
            [ meta, fq2 ? [fq1, fq2] : [fq1] ]
        }
}
