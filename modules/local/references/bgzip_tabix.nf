process BGZIP_TABIX {
    tag "${meta.name}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${vcf}.gz"), path("${vcf}.gz.tbi"), emit: vcf_tbi
    path "versions.yml",                                        emit: versions

    script:
    """
    bgzip -@ ${task.cpus} ${vcf}
    tabix -p vcf ${vcf}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htslib: \$(bgzip --version | head -1 | sed 's/.* //')
    END_VERSIONS
    """
}
