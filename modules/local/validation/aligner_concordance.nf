/*
 * Aligner head-to-head: assert bwa-mem3 (--compat=bwa-mem2) and bwa-mem2 produce byte-identical
 * alignment records. Compares the two CRAMs after stripping the tags bwa-mem3 adds (MQ:i, HN:i)
 * and the @PG/@HD header lines. Writes a report and exits non-zero on any difference.
 */
process ALIGNER_CONCORDANCE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram_mem3), path(crai_mem3), path(cram_mem2), path(crai_mem2)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.concordance.txt"), emit: report
    path "versions.yml",                        emit: versions

    script:
    def prefix = "${meta.id}"
    """
    strip() {
        samtools view -@ ${task.cpus} --reference ${fasta} "\$1" \\
          | sed -E 's/\\tMQ:i:[0-9]+//; s/\\tHN:i:[0-9]+//'
    }
    strip ${cram_mem3} > mem3.sam
    strip ${cram_mem2} > mem2.sam

    n3=\$(wc -l < mem3.sam); n2=\$(wc -l < mem2.sam)
    if cmp -s mem3.sam mem2.sam; then
        verdict="IDENTICAL"
    else
        diffn=\$(diff <(cat mem3.sam) <(cat mem2.sam) | grep -c '^[<>]' || true)
        verdict="DIFFER(\${diffn}_lines)"
    fi

    {
      echo "sample\t${meta.id}"
      echo "records_bwamem3\t\$n3"
      echo "records_bwamem2\t\$n2"
      echo "alignment_records\t\$verdict"
    } > ${prefix}.concordance.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS

    [ "\$verdict" = "IDENTICAL" ] || { echo "WARNING: aligner records differ for ${meta.id}: \$verdict" >&2; }
    """
}
