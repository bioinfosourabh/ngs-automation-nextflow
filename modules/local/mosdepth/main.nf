process MOSDEPTH {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/mosdepth:0.3.10--h4e814b3_1'

    publishDir "${params.outdir}/qc/mosdepth/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(bam), path(bai)
    path intervals

    output:
    tuple val(meta), path("*.global.dist.txt"), emit: global_txt
    tuple val(meta), path("*.summary.txt")    , emit: summary_txt
    tuple val(meta), path("*.region.dist.txt"), emit: regions_txt, optional: true
    tuple val(meta), path("*.regions.bed.gz*"), emit: regions_bed, optional: true
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Targeted data: depth over the capture BED. WGS: depth in fixed 1 kb windows.
    def by_cmd = intervals ? "--by ${intervals}" : "--by 1000"
    """
    mosdepth \\
        --threads ${task.cpus} \\
        --no-per-base \\
        --mapq 20 \\
        ${by_cmd} \\
        ${args} \\
        ${prefix} \\
        ${bam}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mosdepth.global.dist.txt ${prefix}.mosdepth.summary.txt ${prefix}.mosdepth.region.dist.txt
    echo "" | gzip > ${prefix}.regions.bed.gz
    touch ${prefix}.regions.bed.gz.csi

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: 0.3.10
    END_VERSIONS
    """
}
