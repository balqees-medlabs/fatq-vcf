/*
 * Mutect2 --mitochondria-mode over one region:
 *   plain alignment   → non-control region
 *   shifted alignment → control region (shifted coordinates)
 */
process MITO_MUTECT2 {
    tag "${meta.id}.${ref_tag}"
    label 'process_medium'

    input:
    tuple val(meta), val(ref_tag), path(bam), path(bai), path(intervals), path(mt_fasta), path(mt_fai), path(mt_dict)

    output:
    tuple val(meta), val(ref_tag), path("*.vcf.gz"), path("*.vcf.gz.tbi"), path("*.stats"), emit: vcf
    path "versions.yml",                                                                    emit: versions

    script:
    def args   = task.ext.args ?: ''
    def mem    = task.memory.toGiga() - 2
    def prefix = "${meta.id}.chrM.${ref_tag}"
    """
    gatk --java-options "-Xmx${mem}g" Mutect2 \\
        -R ${mt_fasta} -I ${bam} -L ${intervals} ${args} \\
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \\
        --read-filter MateUnmappedAndUnmappedReadFilter \\
        -O ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
