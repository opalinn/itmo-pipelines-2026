process BWAmem {

    publishDir "${params.outdir}/mapped", mode: 'copy'

    input:
        tuple val(sample_id), path(r1), path(r2), path(reference)

    output:
        tuple val(sample_id),
              path("${sample_id}.sorted.bam"),
              path("${sample_id}.sorted.bam.bai"),
              path("${sample_id}.flagstat.txt")

    script:
    """
    bwa index ${reference}
    bwa mem -t ${task.cpus} ${reference} ${r1} ${r2} | \
        samtools view -bS - | \
        samtools sort -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt
    """
}
