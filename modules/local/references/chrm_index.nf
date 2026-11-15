// .fai + .dict + bwa-mem2 index for the 16.5 kb chrM / shifted chrM references (seconds).
process CHRM_INDEX {
    tag "${meta.name}"
    label 'process_low'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta}"), path("${fasta}.fai"), path("${fasta.baseName}.dict"), path("${fasta}.{0123,amb,ann,bwt.2bit.64,pac}"), emit: ref

    script:
    """
    samtools faidx ${fasta}
    samtools dict ${fasta} -o ${fasta.baseName}.dict
    bwa-mem2 index ${fasta}
    """
}
