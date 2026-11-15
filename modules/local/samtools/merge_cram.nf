// Merge recalibrated CRAM shards into one indexed per-sample CRAM (the published analysis-ready file).
process SAMTOOLS_MERGE_CRAM {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(crams)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.recal.cram"), path("*.recal.cram.crai"), emit: cram
    path "versions.yml",                                               emit: versions

    script:
    def prefix = "${meta.id}"
    """
    samtools merge -@ ${task.cpus} --reference ${fasta} --write-index \\
        -o ${prefix}.recal.cram ${crams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
