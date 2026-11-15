process GATK4_GENOTYPEGVCFS {
    tag "chunk_${chunk}"
    label 'process_medium'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    // `input` is either a GenomicsDB workspace directory or a combined GVCF.
    tuple val(chunk), path(input), path(interval)
    path fasta
    path fai
    path dict
    path dbsnp
    path dbsnp_tbi

    output:
    tuple val(chunk), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def input_cmd = input.name.endsWith('.vcf.gz') ? "${input}" : "gendb://${input}"
    def index_cmd = input.name.endsWith('.vcf.gz') ? "gatk IndexFeatureFile --input ${input}" : ''
    def dbsnp_cmd = dbsnp ? "--dbsnp ${dbsnp}" : ''
    """
    ${index_cmd}

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GenotypeGVCFs \\
        --reference ${fasta} \\
        --variant ${input_cmd} \\
        --intervals ${interval} \\
        ${dbsnp_cmd} \\
        --output cohort.${chunk}.vcf.gz \\
        --tmp-dir . \\
        ${args}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > cohort.${chunk}.vcf.gz
    touch cohort.${chunk}.vcf.gz.tbi

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
