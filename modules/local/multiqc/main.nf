process MULTIQC {
    tag 'run'
    label 'process_medium'

    input:
    path  'multiqc_files/*'
    path  config
    path  versions

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data
    path "versions.yml",        emit: versions

    script:
    """
    multiqc --force --config ${config} multiqc_files ${versions}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """
}
