process BWA_INDEX {
    tag "${fasta.name}"
    label 'process_medium'
    container 'quay.io/biocontainers/bwa:0.7.18--he4a0461_0'

    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' || !params.save_reference ? null : fn }

    input:
    path fasta

    output:
    path "bwa"          , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir bwa
    bwa index -p bwa/${fasta.name} ${fasta}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """

    stub:
    """
    mkdir bwa
    for ext in amb ann bwt pac sa; do touch bwa/${fasta.name}.\$ext; done

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.18-r1243-dirty
    END_VERSIONS
    """
}
