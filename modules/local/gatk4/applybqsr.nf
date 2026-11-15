// ApplyBQSR over one interval shard → recalibrated CRAM shard.
process GATK4_APPLYBQSR {
    tag "${meta.id}.${intervals.baseName}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai), path(table), path(intervals)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.recal.cram"), emit: cram
    path "versions.yml",                   emit: versions

    script:
    def args = task.ext.args ?: ''
    def mem  = task.memory.toGiga() - 2
    """
    gatk --java-options "-Xmx${mem}g" ApplyBQSR \\
        -R ${fasta} -I ${cram} -L ${intervals} --bqsr-recal-file ${table} ${args} \\
        -O ${meta.id}.${intervals.baseName}.recal.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
