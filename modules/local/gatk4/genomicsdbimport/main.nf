process GATK4_GENOMICSDBIMPORT {
    tag "chunk_${chunk}"
    label 'process_medium'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(chunk), path(gvcfs), path(tbis), path(interval)

    output:
    tuple val(chunk), path("genomicsdb_${chunk}"), path(interval), emit: genomicsdb
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def inputs    = (gvcfs instanceof List ? gvcfs : [gvcfs]).collect { g -> "--variant ${g}" }.join(' ')
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GenomicsDBImport \\
        ${inputs} \\
        --intervals ${interval} \\
        --merge-input-intervals \\
        --batch-size 50 \\
        --reader-threads ${task.cpus} \\
        --genomicsdb-workspace-path genomicsdb_${chunk} \\
        --tmp-dir . \\
        ${args}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    """
    mkdir genomicsdb_${chunk}
    touch genomicsdb_${chunk}/callset.json

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
