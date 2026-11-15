process GATK4_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'
    container 'broadinstitute/gatk:4.6.1.0'

    publishDir "${params.outdir}/preprocessing/markduplicates/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : (fn.endsWith('.metrics.txt') || params.save_markdup ? fn : null) }

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.md.bam"), path("*.md.bai"), emit: bam
    tuple val(meta), path("*.metrics.txt")              , emit: metrics
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    // One sample may have several lane-level BAMs; MarkDuplicates merges them.
    def bam_list  = bams instanceof List ? bams : [bams]
    def inputs    = bam_list.collect { b -> "--INPUT ${b}" }.join(' ')
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MarkDuplicates \\
        ${inputs} \\
        --OUTPUT ${prefix}.md.bam \\
        --METRICS_FILE ${prefix}.md.metrics.txt \\
        --ASSUME_SORT_ORDER coordinate \\
        --CREATE_INDEX true \\
        --TMP_DIR . \\
        ${args}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.md.bam ${prefix}.md.bai ${prefix}.md.metrics.txt

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
