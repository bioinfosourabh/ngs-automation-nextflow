process GATK4_MERGEMUTECTSTATS {
    tag "${meta.id}"
    label 'process_single'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*.vcf.gz.stats"), emit: stats
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def inputs    = (stats instanceof List ? stats : [stats]).collect { s -> "--stats ${s}" }.join(' ')
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MergeMutectStats \\
        ${inputs} \\
        --output ${prefix}.vcf.gz.stats

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz.stats

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
