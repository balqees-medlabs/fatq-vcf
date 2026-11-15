/*
 * Per sample: merge lane BAMs (if >1), mark duplicates (optical-aware), write indexed CRAM.
 * Input BAMs must be coordinate-sorted with fixmate -m tags (see FASTP_BWAMEM).
 * Stats file (-f) is parsed by MultiQC's samtools markdup module.
 */
process SAMTOOLS_MARKDUP {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(bams)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.cram"), path("*.crai"), emit: cram
    tuple val(meta), path("*.markdup.txt"),          emit: stats
    path "versions.yml",                             emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.id}"
    def input  = bams instanceof List && bams.size() > 1
                    ? "samtools merge -@ ${task.cpus} -u -o - ${bams.join(' ')} |"
                    : "cat ${bams} |"
    """
    ${input} samtools markdup -@ ${task.cpus} ${args} \\
        -T ${prefix}.mdtmp -f ${prefix}.markdup.txt \\
        --reference ${fasta} --output-fmt cram,version=3.0 --write-index \\
        - ${prefix}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
