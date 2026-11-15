// WES only: on-target rate, fold-80, % target bases ≥ 20x/30x. Needs the capture BED as interval_list.
process PICARD_COLLECTHSMETRICS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai), path(bed)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.hs_metrics.txt"), emit: metrics
    path "versions.yml",                      emit: versions

    script:
    def prefix = "${meta.id}"
    def mem    = task.memory.toGiga() - 2
    """
    gatk --java-options "-Xmx${mem}g" BedToIntervalList -I ${bed} -SD ${dict} -O targets.interval_list

    gatk --java-options "-Xmx${mem}g" CollectHsMetrics \\
        -I ${cram} -R ${fasta} \\
        --BAIT_INTERVALS targets.interval_list --TARGET_INTERVALS targets.interval_list \\
        --COVERAGE_CAP 2000 \\
        -O ${prefix}.hs_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
