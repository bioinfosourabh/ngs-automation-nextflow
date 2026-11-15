process BWAMEM2_MEM {
    tag "${meta.id}"
    label 'process_high'
    // bwa-mem2 2.2.1 + samtools 1.19.2 (BioContainers multi-package image)
    container 'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0'

    publishDir "${params.outdir}/preprocessing/aligned/${meta.sample}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' || !params.save_aligned ? null : fn }

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sort_threads = Math.max(1, task.cpus.intdiv(4))
    """
    INDEX=\$(find -L ./ -name "*.amb" | sed 's/\\.amb\$//')

    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -K 100000000 \\
        -Y \\
        -R "${meta.read_group}" \\
        ${args} \\
        \$INDEX \\
        ${reads[0]} ${reads[1]} \\
        | samtools sort -@ ${sort_threads} -T ${prefix}.tmp -o ${prefix}.sorted.bam -

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>/dev/null | tail -n 1)
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.bam

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: 2.2.1
        samtools: 1.19.2
    END_VERSIONS
    """
}
