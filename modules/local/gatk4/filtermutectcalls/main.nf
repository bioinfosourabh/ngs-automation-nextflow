process GATK4_FILTERMUTECTCALLS {
    tag "${meta.id}"
    label 'process_low'
    container 'broadinstitute/gatk:4.6.1.0'

    publishDir "${params.outdir}/variant_calling/somatic/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats), path(artifactprior), path(contamination), path(segmentation)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.filteringStats.tsv")                            , emit: filtering_stats
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}.mutect2"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def orientation_cmd   = artifactprior ? "--ob-priors ${artifactprior}" : ''
    def contamination_cmd = contamination ? "--contamination-table ${contamination}" : ''
    def segmentation_cmd  = segmentation  ? "--tumor-segmentation ${segmentation}" : ''
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        FilterMutectCalls \\
        --reference ${fasta} \\
        --variant ${vcf} \\
        --stats ${stats} \\
        ${orientation_cmd} \\
        ${contamination_cmd} \\
        ${segmentation_cmd} \\
        --filtering-stats ${prefix}.filteringStats.tsv \\
        --output ${prefix}.filtered.vcf.gz \\
        --tmp-dir . \\
        ${args}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.mutect2"
    """
    echo "" | gzip > ${prefix}.filtered.vcf.gz
    touch ${prefix}.filtered.vcf.gz.tbi ${prefix}.filteringStats.tsv

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
