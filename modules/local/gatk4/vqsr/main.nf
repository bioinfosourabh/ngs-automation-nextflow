// Variant quality score recalibration for a joint-called cohort VCF.
// Runs VariantRecalibrator and ApplyVQSR for indels first and then SNPs,
// following the order used in the GATK joint-genotyping reference workflow.
process GATK4_VQSR {
    tag "${meta.id}"
    label 'process_medium'
    container 'broadinstitute/gatk:4.6.1.0'

    publishDir "${params.outdir}/variant_calling/germline/vqsr", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fai
    path dict
    path resources       // hapmap, omni, 1000G, mills, dbsnp (+ optional axiom), with their .tbi files
    val  resource_names  // map: hapmap/omni/g1k/mills/dbsnp/axiom -> file name

    output:
    tuple val(meta), path("*.vqsr.vcf.gz"), path("*.vqsr.vcf.gz.tbi"), emit: vcf
    path "*.{recal,tranches}*"                                       , emit: models
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def snp_args   = task.ext.snp_args   ?: ''
    def indel_args = task.ext.indel_args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def avail_mem  = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def r          = resource_names
    def axiom_cmd  = r.axiom ? "--resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${r.axiom}" : ''
    def tranches   = ['100.0', '99.95', '99.9', '99.8', '99.6', '99.5', '99.4', '99.3', '99.0', '98.0', '97.0', '90.0'].collect { t -> "-tranche ${t}" }.join(' ')
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        VariantRecalibrator \\
        --reference ${fasta} \\
        --variant ${vcf} \\
        --mode INDEL \\
        --resource:mills,known=false,training=true,truth=true,prior=12 ${r.mills} \\
        ${axiom_cmd} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2 ${r.dbsnp} \\
        -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR \\
        ${tranches} \\
        --max-gaussians 4 \\
        --output ${prefix}.indel.recal \\
        --tranches-file ${prefix}.indel.tranches \\
        --tmp-dir . \\
        ${indel_args}

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        VariantRecalibrator \\
        --reference ${fasta} \\
        --variant ${vcf} \\
        --mode SNP \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15 ${r.hapmap} \\
        --resource:omni,known=false,training=true,truth=true,prior=12 ${r.omni} \\
        --resource:1000G,known=false,training=true,truth=false,prior=10 ${r.g1k} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${r.dbsnp} \\
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR \\
        ${tranches} \\
        --max-gaussians 6 \\
        --output ${prefix}.snp.recal \\
        --tranches-file ${prefix}.snp.tranches \\
        --tmp-dir . \\
        ${snp_args}

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyVQSR \\
        --reference ${fasta} \\
        --variant ${vcf} \\
        --mode INDEL \\
        --recal-file ${prefix}.indel.recal \\
        --tranches-file ${prefix}.indel.tranches \\
        --truth-sensitivity-filter-level ${params.vqsr_indel_filter_level} \\
        --create-output-variant-index true \\
        --output ${prefix}.indel_recalibrated.vcf.gz \\
        --tmp-dir .

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyVQSR \\
        --reference ${fasta} \\
        --variant ${prefix}.indel_recalibrated.vcf.gz \\
        --mode SNP \\
        --recal-file ${prefix}.snp.recal \\
        --tranches-file ${prefix}.snp.tranches \\
        --truth-sensitivity-filter-level ${params.vqsr_snp_filter_level} \\
        --create-output-variant-index true \\
        --output ${prefix}.vqsr.vcf.gz \\
        --tmp-dir .

    rm -f ${prefix}.indel_recalibrated.vcf.gz ${prefix}.indel_recalibrated.vcf.gz.tbi

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vqsr.vcf.gz
    touch ${prefix}.vqsr.vcf.gz.tbi ${prefix}.snp.recal ${prefix}.snp.tranches ${prefix}.indel.recal ${prefix}.indel.tranches

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
