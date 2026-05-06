# Homework 2. Running pipelines

## How to run pipeline


1. Run on local data:

```bash
nextflow run main.nf --mode local --input_reads "data/*_{1,2}.fq" --outdir results
```

2. Run with data downloading

```bash
nextflow run main.nf --mode sra --sample_id <sample_id> --outdir results
```

3. Run with reference
```bash
nextflow run main.nf --mode sra --sample_id <sample_id> --reference <reference.fasta> --outdir results
```

4. Run de novo assembly
```bash
nextflow run main.nf --mode sra --sample_id <sample_id> --outdir results
```
