process ENSEMBLVEP_VEP {
    tag "${meta.id}"
    label 'process_medium'
    // The VEP release must match the cache version (params.vep_cache_version).
    container 'quay.io/biocontainers/ensembl-vep:114.2--pl5321h2a3209d_0'

    publishDir "${params.outdir}/annotation/${meta.analysis}/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path cache

    output:
    tuple val(meta), path("*.vep.vcf.gz"), path("*.vep.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.vep_summary.html")                    , emit: report
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${vcf.name.replaceAll(/\.vcf\.gz$/, '')}"
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${prefix}.vep.vcf.gz \\
        --vcf \\
        --compress_output bgzip \\
        --offline \\
        --cache \\
        --dir_cache ${cache} \\
        --species ${params.vep_species} \\
        --assembly ${params.vep_genome} \\
        --cache_version ${params.vep_cache_version} \\
        --fasta ${fasta} \\
        --fork ${task.cpus} \\
        --stats_file ${prefix}.vep_summary.html \\
        ${args}
    tabix -p vcf ${prefix}.vep.vcf.gz

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$(vep --help 2>&1 | sed -n 's/^ *ensembl-vep *: *//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${vcf.name.replaceAll(/\.vcf\.gz$/, '')}"
    """
    echo "" | gzip > ${prefix}.vep.vcf.gz
    touch ${prefix}.vep.vcf.gz.tbi ${prefix}.vep_summary.html

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: 114.2
    END_VERSIONS
    """
}
