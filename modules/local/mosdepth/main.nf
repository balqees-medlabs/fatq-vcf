// Fast coverage: genome-wide (WGS, 500 bp windows) or per capture target (WES).
process MOSDEPTH {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai), path(bed)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.mosdepth.summary.txt"),        emit: summary
    tuple val(meta), path("*.mosdepth.global.dist.txt"),    emit: global_dist
    tuple val(meta), path("*.mosdepth.region.dist.txt"),    emit: region_dist, optional: true
    tuple val(meta), path("*.regions.bed.gz"),              emit: regions,     optional: true
    tuple val(meta), path("*.thresholds.bed.gz"),           emit: thresholds,  optional: true
    path "versions.yml",                                    emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.id}"
    def by     = bed ? "--by ${bed}" : (meta.seq_type == 'wgs' ? '--by 500' : '')
    """
    mosdepth --threads ${task.cpus} --fasta ${fasta} ${by} --thresholds 1,10,20,30 ${args} ${prefix} ${cram}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version | sed 's/mosdepth //')
    END_VERSIONS
    """
}
