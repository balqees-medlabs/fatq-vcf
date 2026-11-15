// Single-sample genotyping of a GVCF → VCF.
process GATK4_GENOTYPEGVCFS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(gvcf), path(gvcf_tbi)
    path  fasta
    path  fai
    path  dict
    tuple path(dbsnp), path(dbsnp_tbi)

    output:
    tuple val(meta), path("*.genotyped.vcf.gz"), path("*.genotyped.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                         emit: versions

    script:
    def args = task.ext.args ?: ''
    def mem  = task.memory.toGiga() - 2
    """
    gatk --java-options "-Xmx${mem}g" GenotypeGVCFs \\
        -R ${fasta} -V ${gvcf} ${args} --dbsnp ${dbsnp} \\
        -O ${meta.id}.genotyped.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
