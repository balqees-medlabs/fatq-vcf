// Drop low-mappability / extreme-GC / low-count intervals before building the gCNV model.
process GATK4_FILTERINTERVALS {
    tag 'pon'
    label 'process_medium'

    input:
    path  intervals
    path  counts        // all cohort counts hdf5
    path  fasta
    path  fai
    path  dict

    output:
    path "filtered.interval_list", emit: intervals
    path "versions.yml",           emit: versions

    script:
    def mem      = task.memory.toGiga() - 2
    def countArg = counts.collect { "-I ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx${mem}g" FilterIntervals \\
        -L ${intervals} ${countArg} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -O filtered.interval_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
