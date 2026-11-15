/*
 * truvari bench: benchmark SV calls (Manta diploidSV) against a GIAB SV truth set within its BED.
 */
process TRUVARI_BENCH {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(query_vcf), path(query_tbi)
    tuple path(truth_vcf), path(truth_tbi), path(truth_bed)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("truvari_${meta.id}/summary.json"), emit: summary
    path "truvari_${meta.id}",                                emit: dir
    path "versions.yml",                                      emit: versions

    script:
    def args = task.ext.args ?: '--passonly'
    """
    truvari bench \\
        -b ${truth_vcf} -c ${query_vcf} \\
        --includebed ${truth_bed} -f ${fasta} ${args} \\
        -o truvari_${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(truvari version 2>&1 | grep -oP '[0-9.]+' | head -1)
    END_VERSIONS
    """
}
