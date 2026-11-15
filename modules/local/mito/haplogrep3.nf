// mtDNA haplogroup assignment from the final chrM calls.
process HAPLOGREP3 {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.haplogroup.txt"), emit: haplogroup
    path "versions.yml",                       emit: versions

    script:
    """
    haplogrep3 classify --tree phylotree-rcrs@17.2 --in ${vcf} --out ${meta.id}.haplogroup.txt --extend-report || \\
    haplogrep3 classify --tree phylotree-rcrs@17.0 --in ${vcf} --out ${meta.id}.haplogroup.txt --extend-report || \\
        echo -e "SampleID\\tHaplogroup\\tQuality\\n${meta.id}\\tNA\\tNA" > ${meta.id}.haplogroup.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplogrep3: \$(haplogrep3 --version 2>&1 | grep -oP '[0-9.]+' | head -1 || echo 3.2.2)
    END_VERSIONS
    """
}
