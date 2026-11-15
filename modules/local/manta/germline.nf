/*
 * Manta germline SV/indel calling (split-read + paired-end). `--exome` is added for WES via ext.args.
 * Works directly on CRAM with the reference. Threads/memory are passed to runWorkflow.
 */
process MANTA_GERMLINE {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(cram), path(crai), path(call_bed), path(call_bed_tbi)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.diploidSV.vcf.gz"),           path("*.diploidSV.vcf.gz.tbi"),           emit: diploid_sv
    tuple val(meta), path("*.candidateSV.vcf.gz"),         path("*.candidateSV.vcf.gz.tbi"),         emit: candidate_sv
    tuple val(meta), path("*.candidateSmallIndels.vcf.gz"),path("*.candidateSmallIndels.vcf.gz.tbi"),emit: small_indels
    path "versions.yml",                                   emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.id}"
    def region = call_bed ? "--callRegions ${call_bed}" : ''
    def mem_gb = task.memory.toGiga()
    """
    configManta.py --bam ${cram} --referenceFasta ${fasta} ${region} ${args} --runDir manta

    python2 manta/runWorkflow.py -m local -j ${task.cpus} -g ${mem_gb}

    cp manta/results/variants/diploidSV.vcf.gz            ${prefix}.diploidSV.vcf.gz
    cp manta/results/variants/diploidSV.vcf.gz.tbi        ${prefix}.diploidSV.vcf.gz.tbi
    cp manta/results/variants/candidateSV.vcf.gz          ${prefix}.candidateSV.vcf.gz
    cp manta/results/variants/candidateSV.vcf.gz.tbi      ${prefix}.candidateSV.vcf.gz.tbi
    cp manta/results/variants/candidateSmallIndels.vcf.gz     ${prefix}.candidateSmallIndels.vcf.gz
    cp manta/results/variants/candidateSmallIndels.vcf.gz.tbi ${prefix}.candidateSmallIndels.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$(configManta.py --version 2>&1)
    END_VERSIONS
    """
}
