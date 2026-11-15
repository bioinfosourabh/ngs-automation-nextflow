process GATK4_CALCULATECONTAMINATION {
    tag "${meta.id}"
    label 'process_single'
    container 'broadinstitute/gatk:4.6.1.0'

    publishDir "${params.outdir}/variant_calling/somatic/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(tumour_table), path(normal_table)

    output:
    tuple val(meta), path("*.contamination.table"), path("*.segmentation.table"), emit: tables
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def avail_mem   = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def matched_cmd = normal_table ? "--matched-normal ${normal_table}" : ''
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CalculateContamination \\
        --input ${tumour_table} \\
        ${matched_cmd} \\
        --tumor-segmentation ${prefix}.segmentation.table \\
        --output ${prefix}.contamination.table

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.contamination.table ${prefix}.segmentation.table

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
