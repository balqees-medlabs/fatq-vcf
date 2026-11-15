/*
 * Liftover the shifted (control-region) calls back to standard chrM coordinates, merge with the
 * non-control calls, merge the Mutect2 stats, then FilterMutectCalls --mitochondria-mode.
 *
 * Input tuple carries both ref variants' VCFs + stats:
 *   plain   = non-control-region calls (standard coords)
 *   shifted = control-region calls (shifted coords, lifted back here)
 */
process MITO_MERGE_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(plain_vcf), path(plain_tbi), path(plain_stats),
                     path(shifted_vcf), path(shifted_tbi), path(shifted_stats)
    path  mt_fasta
    path  mt_fai
    path  mt_dict
    path  shiftback_chain
    path  blacklist

    output:
    tuple val(meta), path("*.chrM.final.vcf.gz"), path("*.chrM.final.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                           emit: versions

    script:
    def args = task.ext.args ?: ''
    def mem  = task.memory.toGiga() - 2
    def pfx  = "${meta.id}"
    """
    # Shifted control-region calls → standard chrM coordinates
    gatk --java-options "-Xmx${mem}g" LiftoverVcf \\
        -I ${shifted_vcf} -O ${pfx}.shifted_back.vcf.gz \\
        -R ${mt_fasta} --CHAIN ${shiftback_chain} \\
        --REJECT ${pfx}.rejected.vcf

    # Merge non-control + lifted control region calls
    gatk --java-options "-Xmx${mem}g" MergeVcfs \\
        -I ${pfx}.shifted_back.vcf.gz -I ${plain_vcf} -O ${pfx}.chrM.merged.vcf.gz

    # Merge Mutect2 stats (required by FilterMutectCalls)
    gatk --java-options "-Xmx${mem}g" MergeMutectStats \\
        --stats ${plain_stats} --stats ${shifted_stats} -O ${pfx}.chrM.merged.stats

    gatk --java-options "-Xmx${mem}g" FilterMutectCalls \\
        -R ${mt_fasta} -V ${pfx}.chrM.merged.vcf.gz ${args} \\
        --stats ${pfx}.chrM.merged.stats \\
        -O ${pfx}.chrM.filtered.vcf.gz

    # Blacklist known artifact sites
    gatk --java-options "-Xmx${mem}g" VariantFiltration \\
        -R ${mt_fasta} -V ${pfx}.chrM.filtered.vcf.gz \\
        --mask ${blacklist} --mask-name "blacklisted_site" \\
        -O ${pfx}.chrM.final.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
    END_VERSIONS
    """
}
