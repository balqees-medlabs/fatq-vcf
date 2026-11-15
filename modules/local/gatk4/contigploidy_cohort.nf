// Build the contig-ploidy model from the PoN cohort (COHORT mode).
process GATK4_CONTIGPLOIDY_COHORT {
    tag 'pon'
    label 'process_high'

    input:
    path  counts
    path  intervals

    output:
    path "ploidy-model", emit: model
    path "ploidy-calls", emit: calls
    path "versions.yml", emit: versions

    script:
    def mem      = task.memory.toGiga() - 4
    def countArg = counts.collect { "-I ${it}" }.join(' ')
    """
    # Default priors table shipped with GATK (contig-ploidy-priors). A kit-specific table can
    # replace this via params in future; the standard human priors work for WES.
    cat > priors.tsv <<'TSV'
CONTIG_NAME	PLOIDY_PRIOR_0	PLOIDY_PRIOR_1	PLOIDY_PRIOR_2	PLOIDY_PRIOR_3
chr1	0.01	0.01	0.97	0.01
chrX	0.01	0.49	0.49	0.01
chrY	0.50	0.50	0.00	0.00
TSV

    gatk --java-options "-Xmx${mem}g" DetermineGermlineContigPloidy \\
        -L ${intervals} ${countArg} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --contig-ploidy-priors priors.tsv \\
        --output . --output-prefix ploidy
    # GATK writes ploidy-model/ and ploidy-calls/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
