process FastQC {

    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(qc_type), val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path("*_fastqc.html")

    script:
    """
    fastqc --threads ${task.cpus} ${r1} ${r2}
    """
}