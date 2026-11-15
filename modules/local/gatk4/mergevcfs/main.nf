process GATK4_MERGEVCFS {
    tag "${meta.id}"
    label 'process_low'
    container 'broadinstitute/gatk:4.6.1.0'

    publishDir "${params.outdir}/variant_calling/${meta.analysis}/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(vcfs)
    path dict
    val suffix      // appended to meta.id to form the output name, e.g. '.raw'

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}${suffix}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def inputs    = (vcfs instanceof List ? vcfs : [vcfs]).collect { v -> "--INPUT ${v}" }.join(' ')
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MergeVcfs \\
        ${inputs} \\
        --SEQUENCE_DICTIONARY ${dict} \\
        --OUTPUT ${prefix}.vcf.gz \\
        --TMP_DIR . \\
        ${args}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}${suffix}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
