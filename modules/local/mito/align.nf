/*
 * Re-align the chrM uBAM to one reference (standard chrM OR shifted-by-8000 chrM).
 * SamToFastq | bwa-mem2 mem | MergeBamAlignment (restore read-group/tags) | MarkDuplicates | sort.
 * `ref_tag` (plain | shifted) labels outputs so the two variants stay distinct downstream.
 */
process MITO_ALIGN {
    tag "${meta.id}.${ref_tag}"
    label 'process_medium'

    input:
    tuple val(meta), path(ubam), val(ref_tag), path(mt_fasta), path(mt_fai), path(mt_dict), path(mt_index)

    output:
    tuple val(meta), val(ref_tag), path("*.bam"), path("*.bai"), emit: bam
    path "versions.yml",                                         emit: versions

    script:
    def mem    = task.memory.toGiga() - 3
    def prefix = "${meta.id}.chrM.${ref_tag}"
    """
    # uBAM → FASTQ (interleaved) → bwa-mem2 → merge back with uBAM metadata
    gatk --java-options "-Xmx2g" SamToFastq -I ${ubam} -F /dev/stdout --INTERLEAVE true --NON_PF true \\
      | bwa-mem2 mem -K 100000000 -Y -t ${task.cpus} -p ${mt_fasta} - 2> >(tee ${prefix}.bwa.err >&2) \\
      | gatk --java-options "-Xmx${mem}g" MergeBamAlignment \\
            --ALIGNED_BAM /dev/stdin --UNMAPPED_BAM ${ubam} --OUTPUT merged.bam \\
            -R ${mt_fasta} --SORT_ORDER unsorted --IS_BISULFITE_SEQUENCE false \\
            --ALIGNED_READS_ONLY false --CLIP_ADAPTERS false --MAX_RECORDS_IN_RAM 2000000 \\
            --ATTRIBUTES_TO_RETAIN X0 --MAX_INSERTIONS_OR_DELETIONS -1 \\
            --PRIMARY_ALIGNMENT_STRATEGY MostDistant --UNMAPPED_READ_STRATEGY COPY_TO_TAG \\
            --ALIGNER_PROPER_PAIR_FLAGS true --UNMAP_CONTAMINANT_READS false --ADD_PG_TAG_TO_READS false

    gatk --java-options "-Xmx${mem}g" MarkDuplicates \\
        -I merged.bam -O md.bam -M ${prefix}.metrics \\
        --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
        --ASSUME_SORT_ORDER queryname --CLEAR_DT false --ADD_PG_TAG_TO_READS false

    gatk --java-options "-Xmx${mem}g" SortSam -I md.bam -O ${prefix}.bam \\
        --SORT_ORDER coordinate --CREATE_INDEX true --MAX_RECORDS_IN_RAM 300000
    mv ${prefix}.bai ${prefix}.bam.bai 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>/dev/null | grep -oP 'v\\K[0-9.]+' | head -1)
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -1)
    END_VERSIONS
    """
}
