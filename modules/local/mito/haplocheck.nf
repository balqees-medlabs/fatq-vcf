// mtDNA contamination estimate from the final chrM calls.
process HAPLOCHECK {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.haplocheck.txt"),     emit: report
    tuple val(meta), path("*.haplocheck.raw.txt"), emit: raw, optional: true
    path "versions.yml",                           emit: versions

    script:
    """
    haplocheck --out ${meta.id}.haplocheck.txt ${vcf} || \\
        echo -e "Sample\\tContamination Status\\tContamination Level\\n${meta.id}\\tNA\\tNA" > ${meta.id}.haplocheck.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplocheck: \$(haplocheck --version 2>&1 | grep -oP '[0-9.]+' | head -1 || echo 1.3.3)
    END_VERSIONS
    """
}
