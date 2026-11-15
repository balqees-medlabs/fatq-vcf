// Cross-sample contamination (FREEMIX) using the 1000G b38 SVD shipped in the container.
process VERIFYBAMID2 {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai)
    path  fasta
    path  fai
    val   svd_prefix

    output:
    tuple val(meta), path("*.selfSM"),    emit: selfsm
    tuple val(meta), path("*.Ancestry"),  emit: ancestry
    path "versions.yml",                  emit: versions

    script:
    def prefix = "${meta.id}"
    """
    verifybamid2 \\
        --SVDPrefix ${svd_prefix} \\
        --Reference ${fasta} \\
        --BamFile ${cram} \\
        --NumThread ${task.cpus} \\
        --Output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid2: \$(verifybamid2 --help 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' | head -1 || echo 2.0.3)
    END_VERSIONS
    """
}
