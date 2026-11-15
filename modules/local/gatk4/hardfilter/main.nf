// GATK generic hard filters, applied separately to SNPs and to indels
// (mixed sites are filtered with the indel thresholds). Failing records are
// kept and marked in the FILTER column.
process GATK4_HARDFILTER {
    tag "${meta.id}"
    label 'process_low'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.hardfiltered.vcf.gz"), path("*.hardfiltered.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def gatk_cmd  = "gatk --java-options \"-Xmx${avail_mem}M -XX:-UsePerfData\""
    """
    ${gatk_cmd} SelectVariants -R ${fasta} -V ${vcf} \\
        --select-type-to-include SNP -O snps.vcf.gz
    ${gatk_cmd} SelectVariants -R ${fasta} -V ${vcf} \\
        --select-type-to-include INDEL --select-type-to-include MIXED -O indels.vcf.gz

    ${gatk_cmd} VariantFiltration -R ${fasta} -V snps.vcf.gz -O snps.filtered.vcf.gz \\
        --filter-name QD2        --filter-expression "QD < 2.0" \\
        --filter-name QUAL30     --filter-expression "QUAL < 30.0" \\
        --filter-name SOR3       --filter-expression "SOR > 3.0" \\
        --filter-name FS60       --filter-expression "FS > 60.0" \\
        --filter-name MQ40       --filter-expression "MQ < 40.0" \\
        --filter-name MQRankSum-12.5     --filter-expression "MQRankSum < -12.5" \\
        --filter-name ReadPosRankSum-8   --filter-expression "ReadPosRankSum < -8.0"

    ${gatk_cmd} VariantFiltration -R ${fasta} -V indels.vcf.gz -O indels.filtered.vcf.gz \\
        --filter-name QD2        --filter-expression "QD < 2.0" \\
        --filter-name QUAL30     --filter-expression "QUAL < 30.0" \\
        --filter-name FS200      --filter-expression "FS > 200.0" \\
        --filter-name ReadPosRankSum-20  --filter-expression "ReadPosRankSum < -20.0"

    ${gatk_cmd} MergeVcfs -I snps.filtered.vcf.gz -I indels.filtered.vcf.gz \\
        -O ${prefix}.hardfiltered.vcf.gz

    rm -f snps.* indels.*

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.hardfiltered.vcf.gz
    touch ${prefix}.hardfiltered.vcf.gz.tbi

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
