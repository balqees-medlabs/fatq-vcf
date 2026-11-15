// somalier extract (per sample) — sites genotypes for sex inference, relatedness and swap checks.
process SOMALIER_EXTRACT {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(cram), path(crai)
    path  fasta
    path  fai
    path  sites

    output:
    tuple val(meta), path("*.somalier"), emit: extract
    path "versions.yml",                 emit: versions

    script:
    """
    somalier extract -d . --sites ${sites} -f ${fasta} ${cram}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(somalier --version 2>&1 | sed 's/somalier version: //')
    END_VERSIONS
    """
}

// somalier relate over every sample in the run — writes samples.tsv (inferred sex) + pairs.tsv.
process SOMALIER_RELATE {
    tag 'all_samples'
    label 'process_low'

    input:
    path extracts

    output:
    path "somalier.samples.tsv", emit: samples
    path "somalier.pairs.tsv",   emit: pairs
    path "somalier.html",        emit: html
    path "versions.yml",         emit: versions

    script:
    """
    somalier relate --output-prefix somalier ${extracts}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(somalier --version 2>&1 | sed 's/somalier version: //')
    END_VERSIONS
    """
}
