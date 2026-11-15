// samtools stats + flagstat + idxstats on the final CRAM (all parsed by MultiQC).
process SAMTOOLS_STATS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.stats"),    emit: stats
    tuple val(meta), path("*.flagstat"), emit: flagstat
    tuple val(meta), path("*.idxstats"), emit: idxstats
    path "versions.yml",                 emit: versions

    script:
    def prefix = "${meta.id}"
    """
    samtools stats    -@ ${task.cpus} --reference ${fasta} ${cram} > ${prefix}.stats
    samtools flagstat -@ ${task.cpus} ${cram} > ${prefix}.flagstat
    samtools idxstats ${cram} > ${prefix}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
