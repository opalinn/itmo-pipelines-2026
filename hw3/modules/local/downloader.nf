process Downloader {
    publishDir "${params.outdir}/raw_reads", mode: 'copy'
    
    input:
        val sample_id

    output:
        tuple val(sample_id), path("${sample_id}_1.fastq"), path("${sample_id}_2.fastq")

    script:
    """
    fastq-dump --split-files ${sample_id}
    """
}