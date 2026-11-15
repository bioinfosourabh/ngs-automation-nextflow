// Keep only records whose FILTER column is PASS.
process BCFTOOLS_VIEW_PASS {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    publishDir "${params.outdir}/variant_calling/${meta.analysis}/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.pass.vcf.gz"), path("*.pass.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${vcf.name.replaceAll(/\.vcf\.gz$/, '')}"
    """
    bcftools view --apply-filters PASS --output-type z --output ${prefix}.pass.vcf.gz ${vcf}
    bcftools index --tbi ${prefix}.pass.vcf.gz

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${vcf.name.replaceAll(/\.vcf\.gz$/, '')}"
    """
    echo "" | gzip > ${prefix}.pass.vcf.gz
    touch ${prefix}.pass.vcf.gz.tbi

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}
