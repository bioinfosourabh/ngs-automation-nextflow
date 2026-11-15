process GATK4_COMBINEGVCFS {
    tag "chunk_${chunk}"
    label 'process_medium'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(chunk), path(gvcfs), path(tbis), path(interval)
    path fasta
    path fai
    path dict

    output:
    tuple val(chunk), path("*.combined.g.vcf.gz"), path(interval), emit: combined
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def inputs    = (gvcfs instanceof List ? gvcfs : [gvcfs]).collect { g -> "--variant ${g}" }.join(' ')
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CombineGVCFs \\
        --reference ${fasta} \\
        ${inputs} \\
        --intervals ${interval} \\
        --output cohort.${chunk}.combined.g.vcf.gz \\
        --tmp-dir . \\
        ${args}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > cohort.${chunk}.combined.g.vcf.gz

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
