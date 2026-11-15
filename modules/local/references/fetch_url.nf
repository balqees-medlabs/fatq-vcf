// Generic resumable download. Runs on the host (needs curl), no container.
process FETCH_URL {
    tag "${meta.subdir}/${meta.name}"
    label 'process_single'
    label 'process_long'
    container null
    maxForks 4

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path(meta.name), emit: file

    script:
    """
    curl -fsSL --retry 5 --retry-delay 10 -C - -o "${meta.name}" "${url}"
    """
}
