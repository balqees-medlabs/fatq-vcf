/*
 * Per-sample calling intervals as a GATK interval_list:
 *   WGS → the Broad wgs_calling_regions list (or params.intervals override)
 *   WES → capture BED converted + padded by params.wes_padding, intersected with the genome
 */
process INTERVALS_PREPARE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(capture_bed)
    path  fasta_dict
    path  wgs_calling_regions
    path  override_intervals   // params.intervals or []

    output:
    tuple val(meta), path("${meta.id}.calling.interval_list"), emit: intervals
    path "versions.yml",                                        emit: versions

    script:
    def mem = task.memory.toGiga() - 1
    if (meta.seq_type == 'wes')
        """
        gatk --java-options "-Xmx${mem}g" BedToIntervalList -I ${capture_bed} -SD ${fasta_dict} -O targets.interval_list
        gatk --java-options "-Xmx${mem}g" IntervalListTools \\
            -I targets.interval_list --PADDING ${params.wes_padding} \\
            --ACTION CONCAT --UNIQUE true --SORT true \\
            -O ${meta.id}.calling.interval_list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
        END_VERSIONS
        """
    else
        """
        cp ${override_intervals ?: wgs_calling_regions} ${meta.id}.calling.interval_list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
        END_VERSIONS
        """
}
