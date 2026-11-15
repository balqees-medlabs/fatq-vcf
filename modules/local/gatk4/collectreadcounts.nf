// Per-sample read counts over the preprocessed gCNV intervals (WES, case mode).
process GATK4_COLLECTREADCOUNTS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai)
    path  fasta
    path  fai
    path  dict
    path  intervals          // preprocessed interval_list from the PoN

    output:
    tuple val(meta), path("*.counts.hdf5"), emit: counts
    path "versions.yml",                    emit: versions

    script:
    def mem = task.memory.toGiga() - 2
    """
    gatk --java-options "-Xmx${mem}g" CollectReadCounts \\
        -R ${fasta} -I ${cram} -L ${intervals} \\
        --interval-merging-rule OVERLAPPING_ONLY --format HDF5 \\
        -O ${meta.id}.counts.hdf5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
