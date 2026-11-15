process GATK4_HAPLOTYPECALLER {
    tag "${meta.id}:${meta.chunk}"
    label 'process_medium'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(meta), path(bam), path(bai), path(interval)
    path fasta
    path fai
    path dict
    path dbsnp
    path dbsnp_tbi

    output:
    tuple val(meta), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), emit: gvcf
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}.${meta.chunk}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def dbsnp_cmd = dbsnp ? "--dbsnp ${dbsnp}" : ''
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        HaplotypeCaller \\
        --input ${bam} \\
        --reference ${fasta} \\
        --intervals ${interval} \\
        --emit-ref-confidence GVCF \\
        --native-pair-hmm-threads ${task.cpus} \\
        ${dbsnp_cmd} \\
        --output ${prefix}.g.vcf.gz \\
        --tmp-dir . \\
        ${args}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chunk}"
    """
    echo "" | gzip > ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
