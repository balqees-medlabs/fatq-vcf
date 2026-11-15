// STR table for HaplotypeCaller --dragen-mode (via CalibrateDragstrModel). Built once per reference.
process GATK4_COMPOSESTRTABLEFILE {
    tag "${fasta.name}"
    label 'process_medium'
    label 'process_long'

    input:
    path fasta
    path fai
    path dict

    output:
    path "GRCh38.str_table.zip", emit: str_table
    path "versions.yml",         emit: versions

    script:
    def mem = task.memory.toGiga() - 2
    """
    gatk --java-options "-Xmx${mem}g" ComposeSTRTableFile -R ${fasta} -O GRCh38.str_table.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
