process GATK4_SPLITINTERVALS {
    tag "scatter_${scatter_count}"
    label 'process_single'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    path fasta
    path fai
    path dict
    path intervals
    val scatter_count

    output:
    path "scatter/*.interval_list", emit: interval_lists
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def interval_cmd = intervals ? "--intervals ${intervals} --interval-padding ${params.interval_padding}" : ''
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SplitIntervals \\
        --reference ${fasta} \\
        ${interval_cmd} \\
        --scatter-count ${scatter_count} \\
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \\
        --output scatter \\
        ${args}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    """
    mkdir scatter
    for i in \$(seq 0 \$(( ${scatter_count} - 1 ))); do
        touch scatter/\$(printf '%04d' \$i)-scattered.interval_list
    done

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
