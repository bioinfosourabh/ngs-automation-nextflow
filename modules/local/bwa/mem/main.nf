process BWA_MEM {
    tag "${meta.id}"
    label 'process_high'
    // bwa 0.7.18 + samtools 1.20 (BioContainers multi-package image)
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0'

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

    bwa mem \\
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
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.bam

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.18-r1243-dirty
        samtools: 1.20
    END_VERSIONS
    """
}
