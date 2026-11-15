# fatq-vcf mito image: GATK4 + bwa-mem2 + samtools for the chrM realignment + Mutect2 pipeline.
# Build:  docker build --network=host -t fatqvcf/mito:0.1.0 containers/mito
FROM mambaorg/micromamba:2.3.3

LABEL org.opencontainers.image.title="fatq-vcf mito" \
      org.opencontainers.image.description="gatk4 4.6.2.0, bwa-mem2 2.3, samtools 1.25 for the GATK mitochondria pipeline" \
      org.opencontainers.image.source="https://github.com/balqees-medlabs/fatq-vcf"

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        gatk4=4.6.2.0 \
        bwa-mem2=2.3 \
        samtools=1.25 \
        procps-ng \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
USER root
