process GATK4_GETPILEUPSUMMARIES {
    tag "${meta.id}"
    label 'process_low'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path dict
    path variants
    path variants_tbi
    path intervals

    output:
    tuple val(meta), path("*.pileups.table"), emit: table
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    // Restrict the common-variant sites to the target regions when they are given.
    def interval_cmd = intervals ? "--intervals ${intervals} --interval-set-rule INTERSECTION" : ''
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GetPileupSummaries \\
        --input ${bam} \\
        --reference ${fasta} \\
        --variant ${variants} \\
        --intervals ${variants} \\
        ${interval_cmd} \\
        --output ${prefix}.pileups.table \\
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
    touch ${prefix}.pileups.table

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
