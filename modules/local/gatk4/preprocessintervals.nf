// Preprocess a WES capture BED into binned, padded gCNV intervals (shared by PoN build + case).
process GATK4_PREPROCESSINTERVALS {
    tag "${bed.name}"
    label 'process_low'

    input:
    path  bed
    path  fasta
    path  fai
    path  dict

    output:
    path "preprocessed.interval_list", emit: intervals
    path "versions.yml",               emit: versions

    script:
    def mem = task.memory.toGiga() - 1
    """
    gatk --java-options "-Xmx${mem}g" PreprocessIntervals \\
        -R ${fasta} -L ${bed} \\
        --bin-length 0 --padding 250 --interval-merging-rule OVERLAPPING_ONLY \\
        -O preprocessed.interval_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
