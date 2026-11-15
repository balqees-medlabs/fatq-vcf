// bcftools stats on the final filtered VCF (parsed by MultiQC).
process BCFTOOLS_STATS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.bcftools_stats.txt"), emit: stats
    path "versions.yml",                           emit: versions

    script:
    """
    bcftools stats -f PASS,. ${vcf} > ${meta.id}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
