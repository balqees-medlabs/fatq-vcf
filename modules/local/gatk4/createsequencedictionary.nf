process GATK4_CREATESEQUENCEDICTIONARY {
    tag "${fasta.name}"
    label 'process_low'

    input:
    path fasta

    output:
    path "*.dict",       emit: dict
    path "versions.yml", emit: versions

    script:
    def mem = task.memory.toGiga() - 1
    """
    gatk --java-options "-Xmx${mem}g" CreateSequenceDictionary -R ${fasta} -O ${fasta.baseName}.dict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
