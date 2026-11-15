// Apply sample-level QC gates (bin/qc_gates.py); never fails the run, flags the sample.
process QC_GATES {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(mosdepth_summary), path(thresholds), path(markdup), path(selfsm)
    path  somalier_samples

    output:
    tuple val(meta), path("*.qc_summary.tsv"), emit: summary
    path "*_qc_mqc.tsv",                       emit: mqc
    path "versions.yml",                       emit: versions

    script:
    def thr = thresholds ? "--thresholds ${thresholds}" : ''
    def som = somalier_samples ? "--somalier ${somalier_samples}" : ''
    """
    qc_gates.py \\
        --sample ${meta.id} --seq-type ${meta.seq_type} --expected-sex ${meta.sex} \\
        --mosdepth ${mosdepth_summary} ${thr} --markdup ${markdup} --selfsm ${selfsm} ${som} \\
        --max-contamination ${params.max_contamination} \\
        --min-mean-cov-wgs ${params.min_mean_cov_wgs} \\
        --min-target-cov-wes ${params.min_target_cov_wes} \\
        --min-target-frac-wes ${params.min_target_frac_wes} \\
        --max-duplication ${params.max_duplication}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
