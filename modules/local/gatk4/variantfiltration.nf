// Hard-filter flags (DRAGEN single-sample QUAL filter by default; see conf/modules.config).
process GATK4_VARIANTFILTRATION {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                       emit: versions

    script:
    def args   = task.ext.args ?: ''
    def mem    = task.memory.toGiga() - 1
    def prefix = task.ext.prefix ?: "${meta.id}.filtered"
    """
    gatk --java-options "-Xmx${mem}g" VariantFiltration \\
        -R ${fasta} -V ${vcf} ${args} \\
        -O ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
