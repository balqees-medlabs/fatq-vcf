process SAMTOOLS_FAIDX {
    tag "${fasta.name}"
    label 'process_single'

    input:
    path fasta

    output:
    path "${fasta}.fai", emit: fai
    path "versions.yml", emit: versions

    script:
    """
    samtools faidx ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
