FROM continuumio/miniconda3

WORKDIR /app

RUN conda install -y -c conda-forge -c bioconda \
    python=3.11 \
    fastqc \
    fastp \
    bwa \
    samtools \
    bcftools \
    spades \
    matplotlib \
    trimmomatic \
    && conda clean -afy

ENV PATH=/opt/conda/bin:$PATH