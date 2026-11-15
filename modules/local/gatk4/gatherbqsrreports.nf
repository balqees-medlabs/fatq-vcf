// Gather scattered BaseRecalibrator tables into one per-sample recalibration model.
process GATK4_GATHERBQSRREPORTS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(tables)

    output:
    tuple val(meta), path("*.recal.table"), emit: table
    path "versions.yml",                    emit: versions

    script:
    def mem   = task.memory.toGiga() - 1
    def input = tables.collect { "-I ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx${mem}g" GatherBQSRReports ${input} -O ${meta.id}.recal.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
