// Build the germline-CNV model from the PoN cohort (COHORT mode). Single scatter for typical WES.
process GATK4_GERMLINECNV_COHORT {
    tag 'pon'
    label 'process_high'
    label 'process_long'

    input:
    path  counts
    path  intervals
    path  ploidy_calls

    output:
    path "cohort-model",  emit: model
    path "cohort-calls",  emit: calls
    path "versions.yml",  emit: versions

    script:
    def mem      = task.memory.toGiga() - 6
    def countArg = counts.collect { "-I ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx${mem}g" GermlineCNVCaller \\
        --run-mode COHORT \\
        -L ${intervals} ${countArg} \\
        --contig-ploidy-calls ${ploidy_calls} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --output . --output-prefix cohort

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
