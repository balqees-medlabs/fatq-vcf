/*
 * AnnotSV: annotate SV/CNV calls (genes, known pathogenic SVs, ACMG class). Optional (needs the
 * AnnotSV annotation bundle at params.annotsv_dir); skipped when not provided.
 */
process ANNOTSV {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(sv_vcf)
    path  annotsv_dir

    output:
    tuple val(meta), path("*.annotated.tsv"), emit: tsv, optional: true
    path "versions.yml",                      emit: versions

    script:
    def prefix = "${meta.id}"
    """
    AnnotSV -SVinputFile ${sv_vcf} -annotationsDir ${annotsv_dir} \\
        -genomeBuild GRCh38 -outputFile ${prefix}.annotated.tsv -outputDir . || \\
        echo "AnnotSV failed (check annotation bundle)" > ${prefix}.annotsv.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(AnnotSV -version 2>&1 | grep -oP '[0-9.]+' | head -1 || echo 3.5.10)
    END_VERSIONS
    """
}
