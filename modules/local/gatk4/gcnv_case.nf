/*
 * GATK germline CNV, case mode (WES): determine contig ploidy against the PoN ploidy model, call
 * CNVs against the PoN cnv model, then postprocess into per-sample genotyped segment VCFs.
 * Requires a panel of normals built per capture kit (see workflows/pon_build.nf).
 */
process GATK4_GERMLINECNV_CASE {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(counts)
    path  ploidy_model       // PoN DetermineGermlineContigPloidy model dir
    path  cnv_model          // PoN GermlineCNVCaller model dir(s)
    path  dict

    output:
    tuple val(meta), path("*.gcnv.vcf.gz"), path("*.gcnv.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                               emit: versions

    script:
    def mem = task.memory.toGiga() - 4
    """
    gatk --java-options "-Xmx${mem}g" DetermineGermlineContigPloidy \\
        --model ${ploidy_model} -I ${counts} \\
        --output ploidy --output-prefix ${meta.id}

    gatk --java-options "-Xmx${mem}g" GermlineCNVCaller \\
        --run-mode CASE --model ${cnv_model} -I ${counts} \\
        --contig-ploidy-calls ploidy/${meta.id}-calls \\
        --output gcnv --output-prefix ${meta.id}

    gatk --java-options "-Xmx${mem}g" PostprocessGermlineCNVCalls \\
        --model-shard-path ${cnv_model} \\
        --calls-shard-path gcnv/${meta.id}-calls \\
        --contig-ploidy-calls ploidy/${meta.id}-calls \\
        --sample-index 0 \\
        --output-genotyped-intervals ${meta.id}.gcnv.intervals.vcf.gz \\
        --output-genotyped-segments ${meta.id}.gcnv.vcf.gz \\
        --sequence-dictionary ${dict}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
