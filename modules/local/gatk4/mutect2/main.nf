process GATK4_MUTECT2 {
    tag "${meta.id}:${meta.chunk}"
    label 'process_medium'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(meta), path(bams), path(bais), path(interval)
    path fasta
    path fai
    path dict
    path germline_resource
    path germline_resource_tbi
    path pon
    path pon_tbi

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.vcf.gz.stats")                , emit: stats
    tuple val(meta), path("*.f1r2.tar.gz")                 , emit: f1r2
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}.${meta.chunk}"
    def avail_mem  = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def inputs     = (bams instanceof List ? bams : [bams]).collect { b -> "--input ${b}" }.join(' ')
    // Read-group sample name of the matched normal; absent in tumour-only mode.
    def normal_cmd = meta.normal_id ? "--normal-sample ${meta.normal_id}" : ''
    def pon_cmd    = pon ? "--panel-of-normals ${pon}" : ''
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        Mutect2 \\
        --reference ${fasta} \\
        ${inputs} \\
        ${normal_cmd} \\
        --germline-resource ${germline_resource} \\
        ${pon_cmd} \\
        --intervals ${interval} \\
        --f1r2-tar-gz ${prefix}.f1r2.tar.gz \\
        --native-pair-hmm-threads ${task.cpus} \\
        --output ${prefix}.vcf.gz \\
        --tmp-dir . \\
        ${args}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chunk}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi ${prefix}.vcf.gz.stats
    echo "" | gzip > ${prefix}.f1r2.tar.gz

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
