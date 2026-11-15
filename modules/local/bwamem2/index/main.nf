process BWAMEM2_INDEX {
    tag "${fasta.name}"
    label 'process_high_memory'
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' || !params.save_reference ? null : fn }

    input:
    path fasta

    output:
    path "bwamem2"      , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir bwamem2
    bwa-mem2 index -p bwamem2/${fasta.name} ${fasta}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>/dev/null | tail -n 1)
    END_VERSIONS
    """

    stub:
    """
    mkdir bwamem2
    for ext in 0123 amb ann bwt.2bit.64 pac; do touch bwamem2/${fasta.name}.\$ext; done

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: 2.2.1
    END_VERSIONS
    """
}
