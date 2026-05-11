process Fastp {

    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz")

    script:
    """
    fastp --in1 ${r1} --in2 ${r2} \
          --out1 ${sample_id}_R1_trimmed.fastq.gz \
          --out2 ${sample_id}_R2_trimmed.fastq.gz \
          --detect_adapter_for_pe \
          --cut_front --cut_tail --cut_window_size 3 --cut_mean_quality 20 \
          --length_required 36 \
          --thread ${task.cpus}
    """
}