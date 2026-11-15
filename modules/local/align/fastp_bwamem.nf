/*
 * Per-lane: fastp (adapter/poly-G trim + QC report) streamed straight into bwa-mem3 / bwa-mem2,
 * then samtools fixmate -m + coordinate sort. No trimmed FASTQ or SAM ever touches disk.
 *
 * Output is a sorted, NOT-yet-duplicate-marked BAM: duplicates are marked per sample after lane
 * merge (MARKDUP), because duplicates must be detected across all lanes of a library.
 *
 *   ext.args  → fastp
 *   ext.args2 → aligner (bwa mem flags)
 */
process FASTP_BWAMEM {
    tag "${meta.id}:${meta.lane}"
    label 'process_align'

    input:
    tuple val(meta), path(reads)
    path  fasta
    path  fai
    path  bwa_index, stageAs: 'bwa_index/*'

    output:
    tuple val(meta), path("*.sorted.bam"),   emit: bam
    tuple val(meta), path("*.fastp.json"),   emit: fastp_json
    tuple val(meta), path("*.fastp.html"),   emit: fastp_html
    path "versions.yml",                     emit: versions

    script:
    def args    = task.ext.args  ?: ''
    def args2   = task.ext.args2 ?: ''
    def prefix  = "${meta.id}_${meta.lane}"
    def rg      = "@RG\\tID:${meta.sample}.${meta.lane}\\tSM:${meta.sample}\\tLB:${meta.sample}\\tPL:ILLUMINA\\tPU:${meta.lane}"
    def idx     = "bwa_index/${fasta.name}"
    // thread split for a 64-cpu task: fastp ~12%, samtools sort 6, aligner the rest
    def t_fastp = Math.max(2, (task.cpus * 0.12) as int)
    def t_sort  = 6
    def t_aln   = Math.max(4, task.cpus - t_fastp - t_sort)
    def mem_per = Math.max(1, ((task.memory.toGiga() * 0.35) / t_sort) as int)
    def mem_sort= "${mem_per}G"
    def aligner = params.aligner == 'bwamem2'
                    ? "bwa-mem2 mem"
                    : "bwa-mem3 mem" + (params.bwamem3_compat ? " --compat=${params.bwamem3_compat}" : '')
    def fq_in   = meta.single_end ? "--in1 ${reads[0]}" : "--in1 ${reads[0]} --in2 ${reads[1]}"
    def interleaved = meta.single_end ? '' : '-p'
    """
    fastp ${fq_in} --stdout --thread ${t_fastp} ${args} \\
          --json ${prefix}.fastp.json --html ${prefix}.fastp.html --report_title "${prefix}" \\
      | ${aligner} -t ${t_aln} ${args2} -R "${rg}" ${interleaved} ${idx} - \\
      | samtools fixmate -@ 2 -m -u - - \\
      | samtools sort -@ ${t_sort} -m ${mem_sort} -T ${prefix}.sorttmp -l 1 -o ${prefix}.sorted.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        aligner: "${aligner}"
        bwa-mem3: \$(bwa-mem3 version 2>&1 | head -1)
        bwa-mem2: \$(micromamba list -n base 2>/dev/null | awk '\$1=="bwa-mem2"{print \$2}')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
