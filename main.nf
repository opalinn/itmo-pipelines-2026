nextflow.enable.dsl = 2

params.mode        = params.mode ?: 'local'
params.input_reads = params.input_reads ?: null
params.sample_id   = params.sample_id ?: null
params.reference   = params.reference ?: null
params.outdir      = params.outdir ?: "results"

process DownloadSRA {
    conda 'bioconda::sra-tools=3.0.3'
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

process FastQCRawData {
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc_raw", mode: 'copy'
    
    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        path("*_fastqc.{zip,html}")

    script:
    """
    fastqc -t ${task.cpus} -q ${r1} ${r2}
    """
}

process Fastp {
    conda 'bioconda::fastp=0.23.4'
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'
    
    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz")

    script:
    """
    fastp -i ${r1} -I ${r2} \\
          -o ${sample_id}_R1_trimmed.fastq.gz -O ${sample_id}_R2_trimmed.fastq.gz \\
          --detect_adapter_for_pe \\
          --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \\
          --length_required 36 \\
          --thread ${task.cpus}
    """
}

process FastQCTrimmed {
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc_trimmed", mode: 'copy'
    
    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        path("*_fastqc.{zip,html}")

    script:
    """
    fastqc -t ${task.cpus} -q ${r1} ${r2}
    """
}

process SPAdes {
    conda 'bioconda::spades=3.15.5'
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
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.19.2'
    publishDir "${params.outdir}/mapped", mode: 'copy'
    
    input:
        tuple val(sample_id), path(r1), path(r2), path(reference)

    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    bwa index ${reference}
    bwa mem -t ${task.cpus} ${reference} ${r1} ${r2} | \
        samtools view -bS - | \
        samtools sort -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}

process CoveragePlot {
    conda 'bioconda::samtools=1.19.2 conda-forge::python=3.10 conda-forge::matplotlib=3.8.0'
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

plt.figure(figsize=(12, 4))
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
        
        raw_reads_ch = DownloadSRA(params.sample_id)
        
    } else {
        if (!params.input_reads)
            error "Specify --input_reads \"data/*_{1,2}.fq\""
        
        reads_pair_ch = Channel.fromFilePairs(params.input_reads, flat: false, checkIfExists: true)
        raw_reads_ch = reads_pair_ch.map { sample_id, reads -> 
            tuple(sample_id, reads[0], reads[1])
        }
    }

    FastQCRawData(raw_reads_ch)
    
    trimmed_reads_ch = Fastp(raw_reads_ch)
    FastQCTrimmed(trimmed_reads_ch)
    
    if (params.reference) {
        ref_ch = Channel.value(file(params.reference))
        BWAmem_in_ch = trimmed_reads_ch.combine(ref_ch).map { sample_id, r1, r2, ref -> 
            tuple(sample_id, r1, r2, ref)
        }
    } else {
        assembly_ch = SPAdes(trimmed_reads_ch)
        BWAmem_in_ch = trimmed_reads_ch.join(assembly_ch).map { sample_id, r1, r2, contigs -> 
            tuple(sample_id, r1, r2, contigs)
        }
    }
    
    bam_ch = BWAmem(BWAmem_in_ch)
    CoveragePlot(bam_ch)
}