// HaplotypeCaller in GVCF + DRAGEN mode over one interval shard.
process GATK4_HAPLOTYPECALLER {
    tag "${meta.id}.${intervals.baseName}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai), path(dragstr_model), path(intervals)
    path  fasta
    path  fai
    path  dict
    tuple path(dbsnp), path(dbsnp_tbi)

    output:
    tuple val(meta), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), emit: gvcf
    path "versions.yml",                                         emit: versions

    script:
    def args    = task.ext.args ?: ''
    def mem     = task.memory.toGiga() - 2
    def dragstr = dragstr_model ? "--dragstr-params-path ${dragstr_model}" : ''
    """
    gatk --java-options "-Xmx${mem}g" HaplotypeCaller \\
        -R ${fasta} -I ${cram} -L ${intervals} ${args} ${dragstr} \\
        --dbsnp ${dbsnp} \\
        -O ${meta.id}.${intervals.baseName}.g.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
