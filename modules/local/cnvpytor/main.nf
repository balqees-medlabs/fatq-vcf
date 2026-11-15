/*
 * CNVpytor read-depth CNV calling (WGS). GC-corrected depth histograms at kb–Mb bins, calls merged
 * into a VCF. Single sample, no panel of normals. BAF import from the sample VCF refines calls.
 */
process CNVPYTOR {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(cram), path(crai), path(vcf), path(vcf_tbi)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.cnvpytor.vcf"),  emit: vcf
    tuple val(meta), path("*.cnvpytor.calls.tsv"), emit: calls
    path "versions.yml",                      emit: versions

    script:
    def prefix = "${meta.id}"
    def bins   = params.cnvpytor_bins ?: '10000 100000'
    def first  = bins.tokenize(' ')[0]
    """
    export OMP_NUM_THREADS=${task.cpus}

    # Read depth (GC correction from the reference FASTA)
    cnvpytor -root ${prefix}.pytor -rd ${cram} -T ${fasta} -j ${task.cpus}

    # Optional BAF from the sample SNV calls
    if [ -s ${vcf} ]; then
        cnvpytor -root ${prefix}.pytor -snp ${vcf} -sample ${meta.id} -nofilter -j ${task.cpus} || true
    fi

    cnvpytor -root ${prefix}.pytor -his ${bins} -j ${task.cpus}
    cnvpytor -root ${prefix}.pytor -partition ${bins} -j ${task.cpus}

    : > ${prefix}.cnvpytor.calls.tsv
    for b in ${bins}; do
        cnvpytor -root ${prefix}.pytor -call \$b >> ${prefix}.cnvpytor.calls.tsv || true
    done

    # Calls → VCF at the finest bin
    printf "set print_filename ${prefix}.cnvpytor.vcf\\nset Q0_range 0 0.5\\nset p_range 0 0.0001\\nprint calls\\nquit\\n" \\
      | cnvpytor -root ${prefix}.pytor -view ${first}
    [ -s ${prefix}.cnvpytor.vcf ] || echo "##fileformat=VCFv4.2" > ${prefix}.cnvpytor.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(cnvpytor --version 2>&1 | grep -oP '[0-9.]+' | head -1)
    END_VERSIONS
    """
}
