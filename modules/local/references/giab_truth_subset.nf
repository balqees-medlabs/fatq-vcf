// GIAB v4.2.1 HG002 truth restricted to chr20, for the regression test's hap.py check.
process GIAB_TRUTH_SUBSET {
    tag 'GIAB_v4.2.1_chr20'
    label 'process_low'

    input:
    val vcf_url
    val bed_url

    output:
    path "HG002_GRCh38_v4.2.1_chr20.vcf.gz*", emit: vcf
    path "HG002_GRCh38_v4.2.1_chr20.bed",     emit: bed

    script:
    """
    bcftools view -r chr20 -Oz -o HG002_GRCh38_v4.2.1_chr20.vcf.gz "${vcf_url}"
    bcftools index -t HG002_GRCh38_v4.2.1_chr20.vcf.gz
    curl -fsSL "${bed_url}" | awk '\$1=="chr20"' > HG002_GRCh38_v4.2.1_chr20.bed
    """
}
