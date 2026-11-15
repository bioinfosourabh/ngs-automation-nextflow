process BCFTOOLS_NORM {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    publishDir "${params.outdir}/variant_calling/${meta.analysis}/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.norm.vcf.gz"), path("*.norm.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${vcf.name.replaceAll(/\.vcf\.gz$/, '')}"
    """
    # Split multi-allelic records and left-align/normalise indels against the reference.
    bcftools norm \\
        --fasta-ref ${fasta} \\
        --multiallelics -any \\
        --check-ref w \\
        --threads ${task.cpus} \\
        ${args} \\
        --output-type z \\
        --output ${prefix}.norm.vcf.gz \\
        ${vcf}
    bcftools index --tbi ${prefix}.norm.vcf.gz

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${vcf.name.replaceAll(/\.vcf\.gz$/, '')}"
    """
    echo "" | gzip > ${prefix}.norm.vcf.gz
    touch ${prefix}.norm.vcf.gz.tbi

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}
