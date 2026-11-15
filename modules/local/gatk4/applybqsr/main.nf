process GATK4_APPLYBQSR {
    tag "${meta.id}"
    label 'process_low'
    container 'broadinstitute/gatk:4.6.1.0'

    publishDir "${params.outdir}/preprocessing/recalibrated/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(bam), path(bai), path(table)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.recal.bam"), path("*.recal.bai"), emit: bam
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyBQSR \\
        --input ${bam} \\
        --reference ${fasta} \\
        --bqsr-recal-file ${table} \\
        --output ${prefix}.recal.bam \\
        --create-output-bam-index true \\
        --tmp-dir . \\
        ${args}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.recal.bam ${prefix}.recal.bai

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
