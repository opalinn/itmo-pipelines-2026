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