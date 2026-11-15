process BCFTOOLS_STATS {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    publishDir "${params.outdir}/qc/bcftools_stats", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta

    output:
    tuple val(meta), path("*.bcftools_stats.txt"), emit: stats
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${vcf.name.replaceAll(/\.vcf\.gz$/, '')}"
    """
    bcftools stats --fasta-ref ${fasta} ${vcf} > ${prefix}.bcftools_stats.txt

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${vcf.name.replaceAll(/\.vcf\.gz$/, '')}"
    """
    touch ${prefix}.bcftools_stats.txt

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}
