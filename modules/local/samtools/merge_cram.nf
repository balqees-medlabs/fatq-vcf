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
    # CRAM 3.0 (not samtools-1.25 default 3.1): GATK 4.6 htsjdk cannot read CRAM 3.1.
    samtools merge -@ ${task.cpus} --reference ${fasta} \\
        --output-fmt cram,version=3.0 --write-index \\
        -o ${prefix}.recal.cram ${crams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
