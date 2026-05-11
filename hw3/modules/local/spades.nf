process SPAdes {

    publishDir "${params.outdir}/assembly", mode: 'copy'
    
    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path("${sample_id}_contigs.fasta")

    script:
    """
    spades.py -t ${task.cpus} -1 ${r1} -2 ${r2} -o spades_out
    cp spades_out/contigs.fasta ${sample_id}_contigs.fasta
    """
}