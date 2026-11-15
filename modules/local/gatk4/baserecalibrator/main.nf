process GATK4_BASERECALIBRATOR {
    tag "${meta.id}"
    label 'process_low'
    container 'broadinstitute/gatk:4.6.1.0'

    publishDir "${params.outdir}/preprocessing/bqsr/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path dict
    path known_sites
    path known_sites_tbi
    path intervals

    output:
    tuple val(meta), path("*.recal.table"), emit: table
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def sites     = (known_sites instanceof List ? known_sites : [known_sites]).collect { v -> "--known-sites ${v}" }.join(' ')
    def interval_cmd = intervals ? "--intervals ${intervals} --interval-padding ${params.interval_padding}" : ''
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        BaseRecalibrator \\
        --input ${bam} \\
        --reference ${fasta} \\
        ${sites} \\
        ${interval_cmd} \\
        --output ${prefix}.recal.table \\
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
    touch ${prefix}.recal.table

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
