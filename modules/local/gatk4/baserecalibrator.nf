// BaseRecalibrator over one interval shard. Scattered; reports gathered downstream.
process GATK4_BASERECALIBRATOR {
    tag "${meta.id}.${intervals.baseName}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai), path(intervals)
    path  fasta
    path  fai
    path  dict
    tuple path(dbsnp), path(dbsnp_tbi)
    tuple path(mills), path(mills_tbi)
    tuple path(known), path(known_tbi)

    output:
    tuple val(meta), path("*.recal.table"), emit: table
    path "versions.yml",                    emit: versions

    script:
    def args = task.ext.args ?: ''
    def mem  = task.memory.toGiga() - 2
    """
    gatk --java-options "-Xmx${mem}g" BaseRecalibrator \\
        -R ${fasta} -I ${cram} -L ${intervals} ${args} \\
        --known-sites ${dbsnp} --known-sites ${mills} --known-sites ${known} \\
        -O ${meta.id}.${intervals.baseName}.recal.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
