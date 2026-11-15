#!/usr/bin/env nextflow
/*
 * ngs-automation-nextflow
 * Germline and somatic short-variant calling from paired-end Illumina FASTQ
 * following the GATK4 best-practice workflows.
 *
 * https://github.com/bioinfosourabh/ngs-automation-nextflow
 */

include { NGS_VARIANTS                                         } from './workflows/ngs_variants'
include { validateParams; pipelineSummary; completionMessage } from './subworkflows/local/utils_pipeline'

workflow {
    main:
    validateParams("${projectDir}/nextflow_schema.json")
    log.info(pipelineSummary())

    NGS_VARIANTS()

    workflow.onComplete = {
        log.info(completionMessage())
    }
}
