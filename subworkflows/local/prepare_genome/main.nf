//
// Reference preparation: FASTA index, sequence dictionary, aligner index and
// interval scattering. Indexes given on the command line are used as-is;
// missing ones are built.
//
include { SAMTOOLS_FAIDX                 } from '../../../modules/local/samtools/faidx'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/local/gatk4/createsequencedictionary'
include { BWAMEM2_INDEX                  } from '../../../modules/local/bwamem2/index'
include { BWA_INDEX                      } from '../../../modules/local/bwa/index'
include { GATK4_SPLITINTERVALS           } from '../../../modules/local/gatk4/splitintervals'

workflow PREPARE_GENOME {
    main:
    ch_versions = channel.empty()
    ch_fasta    = channel.value(file(params.fasta, checkIfExists: true))

    ch_fai = channel.empty()
    if (params.fasta_fai) {
        ch_fai = channel.value(file(params.fasta_fai, checkIfExists: true))
    } else {
        SAMTOOLS_FAIDX(ch_fasta)
        ch_fai      = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    ch_dict = channel.empty()
    if (params.dict) {
        ch_dict = channel.value(file(params.dict, checkIfExists: true))
    } else {
        GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
        ch_dict     = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }

    ch_index = channel.empty()
    if (params.bwa_index) {
        ch_index = channel.value(file(params.bwa_index, checkIfExists: true))
    } else if (params.aligner == 'bwa-mem2') {
        BWAMEM2_INDEX(ch_fasta)
        ch_index    = BWAMEM2_INDEX.out.index
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    } else {
        BWA_INDEX(ch_fasta)
        ch_index    = BWA_INDEX.out.index
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    }

    // Target regions (BED). An empty list stands for "whole genome".
    ch_intervals = params.intervals ? channel.value(file(params.intervals, checkIfExists: true)) : channel.value([])

    // Scatter the calling regions into chunks: [chunk_id, interval_list]
    GATK4_SPLITINTERVALS(ch_fasta, ch_fai, ch_dict, ch_intervals, params.scatter_count)
    ch_scatter = GATK4_SPLITINTERVALS.out.interval_lists
        .flatten()
        .map { f -> [f.name.tokenize('-')[0], f] }
    ch_versions = ch_versions.mix(GATK4_SPLITINTERVALS.out.versions)

    emit:
    fasta     = ch_fasta
    fai       = ch_fai
    dict      = ch_dict
    index     = ch_index
    intervals = ch_intervals
    scatter   = ch_scatter
    versions  = ch_versions
}
