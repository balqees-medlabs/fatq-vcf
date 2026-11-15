// Gather scattered (G)VCF shards into one sorted, indexed VCF. ext.prefix sets the output name/suffix.
process GATK4_MERGEVCFS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(vcfs), path(tbis)
    path  dict

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                     emit: versions

    script:
    def mem    = task.memory.toGiga() - 1
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = vcfs.collect { "-I ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx${mem}g" MergeVcfs ${input} -D ${dict} -O ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
