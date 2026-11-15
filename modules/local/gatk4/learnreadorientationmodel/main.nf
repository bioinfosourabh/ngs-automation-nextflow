process GATK4_LEARNREADORIENTATIONMODEL {
    tag "${meta.id}"
    label 'process_low'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(meta), path(f1r2)

    output:
    tuple val(meta), path("*.artifact-prior.tar.gz"), emit: artifactprior
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def inputs    = (f1r2 instanceof List ? f1r2 : [f1r2]).collect { f -> "--input ${f}" }.join(' ')
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        LearnReadOrientationModel \\
        ${inputs} \\
        --output ${prefix}.artifact-prior.tar.gz \\
        --tmp-dir .

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.artifact-prior.tar.gz

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
