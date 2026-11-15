/*
 * Extract chrM reads from the analysis CRAM and revert to an unmapped, queryname-sorted BAM
 * (strips old alignment/tags so the reads can be re-aligned to the chrM and shifted-chrM refs).
 */
process MITO_SUBSET_REVERT {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.chrM.unmapped.bam"), emit: ubam
    path "versions.yml",                          emit: versions

    script:
    def mem = task.memory.toGiga() - 3
    """
    # chrM reads → BAM
    gatk --java-options "-Xmx${mem}g" PrintReads \\
        -R ${fasta} -I ${cram} -L chrM \\
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \\
        --read-filter MateUnmappedAndUnmappedReadFilter \\
        -O ${meta.id}.chrM.bam

    # Revert alignment → unmapped, queryname-sorted uBAM
    gatk --java-options "-Xmx${mem}g" RevertSam \\
        -I ${meta.id}.chrM.bam \\
        --SANITIZE true --ATTRIBUTE_TO_CLEAR FT --ATTRIBUTE_TO_CLEAR CO \\
        --RESTORE_ORIGINAL_QUALITIES false --OUTPUT_BY_READGROUP false \\
        --SORT_ORDER queryname \\
        -O ${meta.id}.chrM.unmapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
