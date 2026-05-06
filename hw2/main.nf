nextflow.enable.dsl = 2

params.mode        = params.mode ?: 'local'
params.input_reads = params.input_reads ?: null
params.sample_id   = params.sample_id ?: null
params.reference   = params.reference ?: null
params.outdir      = params.outdir ?: "results"

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

process CoveragePlot {

    publishDir "${params.outdir}/coverage", mode: 'copy'
    
    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}_coverage.png"), path("${sample_id}_depth.txt")

    script:
    """
    samtools depth ${bam} > ${sample_id}_depth.txt
    
    python3 -c "
import matplotlib.pyplot as plt

sample_id = '${sample_id}'
depths = []

with open(f'{sample_id}_depth.txt') as f:
    for line in f:
        depths.append(int(line.split()[2]))

plt.figure(dpi=300)
plt.plot(depths, color='darkred', linewidth=0.5)
plt.title(f'Coverage depth: {sample_id}')
plt.xlabel('Position')
plt.ylabel('Depth')
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(f'{sample_id}_coverage.png', dpi=300)
print(f'Total positions: {len(depths)}')
print(f'Mean depth: {sum(depths)/len(depths):.2f}')
print(f'Min depth: {min(depths)}')
print(f'Max depth: {max(depths)}')
"
    """
}

workflow {
    
    if (params.mode == 'sra') {
        if (!params.sample_id)
            error "Specify --sample_id SRR"
        
        raw_reads = Downloader(params.sample_id)
        
    } else {
        if (!params.input_reads)
            error "Specify --input_reads \"data/*_{1,2}.fq\""
        
        pe_reads = channel.fromFilePairs(params.input_reads, flat: false, checkIfExists: true)
        raw_reads = pe_reads.map { sample_id, reads -> 
            tuple(sample_id, reads[0], reads[1])
        }
    }

    qc_raw = raw_reads.map { sample_id, r1, r2 ->
    tuple('qc_raw', sample_id, r1, r2)
    }
    
    trimmed_reads = Fastp(raw_reads)

    qc_trim = trimmed_reads.map { sample_id, r1, r2 ->
    tuple('qc_trimmed', sample_id, r1, r2)
    }

    FastQC(qc_raw.mix(qc_trim))
    
    if (params.reference) {
        ref_ch = channel.value(file(params.reference))
        mapping_data = trimmed_reads.combine(ref_ch).map { sample_id, r1, r2, ref -> 
            tuple(sample_id, r1, r2, ref)
        }
    } else {
        assembly = SPAdes(trimmed_reads)
        mapping_data = trimmed_reads.join(assembly).map { sample_id, r1, r2, contigs -> 
            tuple(sample_id, r1, r2, contigs)
        }
    }
    
    mapping_result = BWAmem(mapping_data)
    coverage_input = mapping_result.map { sample_id, bam, bai, flagstat ->
        tuple(sample_id, bam)
    }

    CoveragePlot(coverage_input)
}