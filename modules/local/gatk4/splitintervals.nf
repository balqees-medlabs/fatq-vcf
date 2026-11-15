/*
 * Split a sample's calling intervals into N balanced shards for scatter.
 * BALANCING_WITHOUT_INTERVAL_SUBDIVISION keeps whole regions together (needed for assembly callers).
 */
process GATK4_SPLITINTERVALS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(intervals)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("shards/*.interval_list"), emit: shards
    path "versions.yml",                             emit: versions

    script:
    def mem = task.memory.toGiga() - 1
    """
    mkdir shards
    gatk --java-options "-Xmx${mem}g" SplitIntervals \\
        -R ${fasta} -L ${intervals} \\
        --scatter-count ${params.scatter_count} \\
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \\
        -O shards

    # Rename 0000-scattered.interval_list → <id>.0000.interval_list for stable tags
    i=0
    for f in shards/*-scattered.interval_list; do
        printf -v n "%04d" \$i
        mv "\$f" "shards/${meta.id}.\${n}.interval_list"
        i=\$((i+1))
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
