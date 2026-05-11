nextflow.enable.dsl = 2

params.mode        = params.mode ?: null
params.input_reads = params.input_reads ?: null
params.sample_id   = params.sample_id ?: null
params.reference   = params.reference ?: null
params.outdir      = params.outdir ?: "results"

//local modules
include { Downloader } from './modules/local/downloader.nf'
include { FastQC } from './modules/local/fastqc.nf'
include { Fastp } from './modules/local/fastp.nf'
include { SPAdes } from './modules/local/spades.nf'
include { BWAmem } from './modules/local/mapping.nf'
include { CoveragePlot } from './modules/local/coverage.nf'

//nf-core modules
include { SAMTOOLS_FAIDX as SAMToolsIdx} from './modules/nf-core/samtools/faidx/main'
workflow {
    if (params.mode == 'sra') {
        if (!params.sample_id)
            error "Specify --sample_id SRR"
        raw_reads = Downloader(params.sample_id)

    } else {
        if (!params.input_reads)
            error "Specify --input_reads \"data/*_{1,2}.fq\""
        pe_reads = channel.fromFilePairs(
            params.input_reads,
            flat: false,
            checkIfExists: true
        )
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
        reference_ch = channel.value(file(params.reference))
        mapping_data = trimmed_reads.combine(reference_ch).map {
            sample_id, r1, r2, ref ->
            tuple(sample_id, r1, r2, ref)
        }

        reference_for_index = channel.value(
            tuple(
                [id: 'reference'],
                file(params.reference),
                []
            )
        )

    } else {
        assembly = SPAdes(trimmed_reads)
        mapping_data = trimmed_reads.join(assembly).map {
            sample_id, r1, r2, contigs ->
            tuple(sample_id, r1, r2, contigs)
        }

        reference_for_index = assembly.map { sample_id, contigs ->
            tuple(
                [id: sample_id],
                contigs,
                []
            )
        }
    }

    SAMToolsIdx(reference_for_index, false)

    reference_for_index.map { meta, reference, _ignore ->
            tuple(meta, reference)
        }
        .join(SAMToolsIdx.out.fai)
        .map { meta, reference, fai ->
            tuple(meta, reference, fai)
        }

    mapping_result = BWAmem(mapping_data)

    coverage_input = mapping_result.map {
        sample_id, bam, _bai, _flagstat ->

        tuple(sample_id, bam)
    }

    CoveragePlot(coverage_input)

}