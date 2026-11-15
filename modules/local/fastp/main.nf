process FASTP {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    publishDir "${params.outdir}/preprocessing/fastp/${meta.sample}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : (fn.endsWith('.fastq.gz') && !params.save_trimmed ? null : fn) }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    tuple val(meta), path("*.fastp.json")   , emit: json
    tuple val(meta), path("*.fastp.html")   , emit: html
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${prefix}_1.trim.fastq.gz \\
        --out2 ${prefix}_2.trim.fastq.gz \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --detect_adapter_for_pe \\
        --thread ${task.cpus} \\
        ${args}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_1.trim.fastq.gz
    echo "" | gzip > ${prefix}_2.trim.fastq.gz
    echo '{}' > ${prefix}.fastp.json
    touch ${prefix}.fastp.html

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        fastp: 0.23.4
    END_VERSIONS
    """
}
