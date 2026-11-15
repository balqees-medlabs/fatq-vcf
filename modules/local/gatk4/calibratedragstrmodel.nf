// Per-sample DRAGstr model for HaplotypeCaller --dragen-mode (STR-aware indel calling).
process GATK4_CALIBRATEDRAGSTRMODEL {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai)
    path  fasta
    path  fai
    path  dict
    path  str_table

    output:
    tuple val(meta), path("*.dragstr_model.txt"), emit: model
    path "versions.yml",                          emit: versions

    script:
    def mem = task.memory.toGiga() - 2
    """
    gatk --java-options "-Xmx${mem}g" CalibrateDragstrModel \\
        -R ${fasta} -I ${cram} -str ${str_table} \\
        -O ${meta.id}.dragstr_model.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
