/*
 * Build a small, realistic WGS test set: HG002 reads on chr20 + chrM streamed over HTTPS from the
 * indexed GIAB GRCh38 BAM, downsampled to ~30x on chr20, converted back to paired FASTQ.
 * Also writes the matching interval_list and samplesheet for -profile test.
 */
process GIAB_TEST_SUBSET {
    tag 'HG002_chr20_chrM'
    label 'process_medium'
    label 'process_long'

    input:
    val bam_url
    path fasta
    path fai
    path dict

    output:
    path "HG002_chr20_chrM_R{1,2}.fastq.gz", emit: fastq
    path "chr20_chrM.interval_list",          emit: interval_list
    path "samplesheet_test_wgs.csv",          emit: samplesheet
    path "subset_stats.txt",                  emit: stats

    script:
    def target_depth = 30
    def chr20_len    = 64444167
    """
    # 1. Stream only chr20 + chrM records (remote .bai fetched automatically by htslib)
    samtools view -@ ${task.cpus} -b -X -o region.bam "${bam_url}" "${bam_url}.bai" chr20 chrM

    # 2. Measure chr20 depth, derive downsample fraction (guard against zero / SIGPIPE)
    reads=\$(samtools view -c -F 0x904 region.bam chr20)
    readlen=\$(samtools view region.bam chr20 | awk 'NR<=1000{s+=length(\$10); n++} NR>1000{exit} END{print (n>0)? int(s/n) : 150}')
    depth=\$(awk -v r=\$reads -v l=\$readlen -v L=${chr20_len} 'BEGIN{printf "%.1f", (r*l)/L}')
    frac=\$(awk -v d=\$depth -v t=${target_depth} 'BEGIN{f=(d>0)? t/d : 1; if (f>1) f=1; printf "%.4f", f}')
    echo -e "chr20_reads\t\$reads\nchr20_mean_read_len\t\$readlen\nchr20_depth_before\t\$depth\ndownsample_fraction\t\$frac" > subset_stats.txt

    # 3. Downsample (seed 42), restore name order, write paired FASTQ (drop singletons/supplementary)
    samtools view -@ ${task.cpus} -b -s 42\${frac#0} -F 0x900 region.bam \\
      | samtools collate -@ ${task.cpus} -u -O - tmp_collate \\
      | samtools fastq -@ ${task.cpus} -1 HG002_chr20_chrM_R1.fastq.gz -2 HG002_chr20_chrM_R2.fastq.gz \\
                        -0 /dev/null -s /dev/null -n -

    # 4. Interval list + samplesheet
    grep -E '^@' ${dict} > chr20_chrM.interval_list
    printf "chr20\t1\t${chr20_len}\t+\tchr20\nchrM\t1\t16569\t+\tchrM\n" >> chr20_chrM.interval_list

    echo "sample,lane,fastq_1,fastq_2,seq_type,capture_bed,sex" > samplesheet_test_wgs.csv
    echo "HG002_chr20,L001,${projectDir}/tests/data/HG002_chr20_chrM_R1.fastq.gz,${projectDir}/tests/data/HG002_chr20_chrM_R2.fastq.gz,wgs,,male" >> samplesheet_test_wgs.csv
    """
}
