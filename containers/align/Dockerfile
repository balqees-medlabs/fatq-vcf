# fatq-vcf alignment image: fastp + bwa-mem3 + bwa-mem2 (fallback) + samtools, all from Bioconda, pinned.
# fastp streams trimmed reads directly into the aligner, so they share one image.
# Build:  docker build -t fatqvcf/align:0.1.0 containers/align
# Versions must match conf/containers.config and docs/PIPELINE_PLAN.md.
FROM mambaorg/micromamba:2.3.3

LABEL org.opencontainers.image.title="fatq-vcf align" \
      org.opencontainers.image.description="fastp 1.4.0, bwa-mem3 0.17.0, bwa-mem2 2.3, samtools 1.25" \
      org.opencontainers.image.source="https://github.com/balqees-medlabs/fatq-vcf"

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        fastp=1.4.0 \
        bwa-mem3=0.17.0 \
        bwa-mem2=2.3 \
        samtools=1.25 \
        procps-ng \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
# Nextflow runs containers with -u $(id -u); make sure no step needs the micromamba user.
USER root
