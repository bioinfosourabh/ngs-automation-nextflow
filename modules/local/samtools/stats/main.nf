process SAMTOOLS_STATS {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    publishDir "${params.outdir}/qc/samtools/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("*.stats")   , emit: stats
    tuple val(meta), path("*.flagstat"), emit: flagstat
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools stats --threads ${task.cpus} --reference ${fasta} ${bam} > ${prefix}.stats
    samtools flagstat --threads ${task.cpus} ${bam} > ${prefix}.flagstat

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats ${prefix}.flagstat

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
