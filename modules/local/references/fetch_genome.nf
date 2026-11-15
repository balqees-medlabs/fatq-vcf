// Download + decompress the GRCh38 analysis-set FASTA. Runs on the host (needs curl), no container.
process FETCH_GENOME {
    tag 'GRCh38_no_alt_plus_hs38d1'
    label 'process_low'
    label 'process_long'
    container null

    input:
    val url

    output:
    path "*.fna", emit: fasta

    script:
    def gz = url.tokenize('/')[-1]
    """
    curl -fsSL --retry 5 --retry-delay 10 -o ${gz} "${url}"
    gunzip ${gz}
    """
}
