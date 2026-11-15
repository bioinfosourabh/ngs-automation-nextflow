process SAMTOOLS_FAIDX {
    tag "${fasta.name}"
    label 'process_single'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' || !params.save_reference ? null : fn }

    input:
    path fasta

    output:
    path "*.fai"        , emit: fai
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools faidx ${fasta}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.fai

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
