/*
 * bwa-mem2 index — read as-is by bwa-mem3 too (.bwt.2bit.64 + .pac); the .0123 file keeps
 * the bwa-mem2 fallback working. Human genome: ~90 GB RAM, ~1-2 h, ~20 GB on disk.
 */
process BWAMEM2_INDEX {
    tag "${fasta.name}"
    label 'process_high_memory'
    label 'process_long'
    cpus 4

    input:
    path fasta

    output:
    path "${fasta}.{0123,amb,ann,bwt.2bit.64,pac}", emit: index
    path "versions.yml",                             emit: versions

    script:
    """
    bwa-mem2 index ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -1)
    END_VERSIONS
    """
}
