//
// Germline short variants: HaplotypeCaller GVCFs per sample and chunk,
// joint genotyping per chunk (GenomicsDBImport or CombineGVCFs, then
// GenotypeGVCFs), gather, filtering (VQSR or hard filters) and normalisation.
//
include { GATK4_HAPLOTYPECALLER                            } from '../../../modules/local/gatk4/haplotypecaller'
include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_GVCF          } from '../../../modules/local/gatk4/mergevcfs'
include { GATK4_GENOMICSDBIMPORT                           } from '../../../modules/local/gatk4/genomicsdbimport'
include { GATK4_COMBINEGVCFS                               } from '../../../modules/local/gatk4/combinegvcfs'
include { GATK4_GENOTYPEGVCFS                              } from '../../../modules/local/gatk4/genotypegvcfs'
include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_COHORT        } from '../../../modules/local/gatk4/mergevcfs'
include { GATK4_VQSR                                       } from '../../../modules/local/gatk4/vqsr'
include { GATK4_HARDFILTER                                 } from '../../../modules/local/gatk4/hardfilter'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_GERMLINE          } from '../../../modules/local/bcftools/norm'
include { BCFTOOLS_VIEW_PASS as BCFTOOLS_VIEW_PASS_GERMLINE } from '../../../modules/local/bcftools/view'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_GERMLINE        } from '../../../modules/local/bcftools/stats'

workflow GERMLINE_JOINT {
    take:
    ch_bam          // [meta, bam, bai] normal samples
    ch_scatter      // [chunk, interval_list]
    ch_fasta
    ch_fai
    ch_dict
    ch_dbsnp
    ch_dbsnp_tbi
    ch_vqsr_files   // list of VQSR resource VCFs + indexes (empty list for hard filtering)
    vqsr_names      // map of resource name -> file name

    main:
    ch_versions = channel.empty()

    ch_hc_in = ch_bam
        .combine(ch_scatter)
        .map { meta, bam, bai, chunk, interval -> [meta + [chunk: chunk], bam, bai, interval] }
    GATK4_HAPLOTYPECALLER(ch_hc_in, ch_fasta, ch_fai, ch_dict, ch_dbsnp, ch_dbsnp_tbi)
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    // Optional: gather per-sample GVCFs so they can be reused for later joint calls.
    if (params.save_gvcf) {
        ch_sample_gvcf = GATK4_HAPLOTYPECALLER.out.gvcf
            .map { meta, gvcf, _tbi -> [meta.findAll { k, _v -> k != 'chunk' } + [analysis: 'germline'], gvcf] }
            .groupTuple()
            .map { meta, gvcfs -> [meta, gvcfs.sort { g -> g.name }] }
        GATK4_MERGEVCFS_GVCF(ch_sample_gvcf, ch_dict, '.g')
    }

    // Group all samples' GVCFs by chunk: [chunk, [gvcf...], [tbi...], interval_list]
    ch_joint_in = GATK4_HAPLOTYPECALLER.out.gvcf
        .map { meta, gvcf, tbi -> [meta.chunk, gvcf, tbi] }
        .groupTuple()
        .map { chunk, gvcfs, tbis -> [chunk, gvcfs.sort { g -> g.name }, tbis.sort { t -> t.name }] }
        .join(ch_scatter)

    ch_genotype_in = channel.empty()
    if (params.joint_caller == 'genomicsdb') {
        GATK4_GENOMICSDBIMPORT(ch_joint_in)
        ch_genotype_in = GATK4_GENOMICSDBIMPORT.out.genomicsdb
        ch_versions    = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())
    } else {
        GATK4_COMBINEGVCFS(ch_joint_in, ch_fasta, ch_fai, ch_dict)
        ch_genotype_in = GATK4_COMBINEGVCFS.out.combined
        ch_versions    = ch_versions.mix(GATK4_COMBINEGVCFS.out.versions.first())
    }

    GATK4_GENOTYPEGVCFS(ch_genotype_in, ch_fasta, ch_fai, ch_dict, ch_dbsnp, ch_dbsnp_tbi)
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions.first())

    ch_cohort_parts = GATK4_GENOTYPEGVCFS.out.vcf
        .map { _chunk, vcf, _tbi -> vcf }
        .collect()
        .map { vcfs -> [[id: 'joint_germline', analysis: 'germline'], vcfs.sort { v -> v.name }] }
    GATK4_MERGEVCFS_COHORT(ch_cohort_parts, ch_dict, '.raw')
    ch_versions = ch_versions.mix(GATK4_MERGEVCFS_COHORT.out.versions)

    ch_filtered = channel.empty()
    if (params.germline_filter == 'vqsr') {
        GATK4_VQSR(GATK4_MERGEVCFS_COHORT.out.vcf, ch_fasta, ch_fai, ch_dict, ch_vqsr_files, vqsr_names)
        ch_filtered = GATK4_VQSR.out.vcf
        ch_versions = ch_versions.mix(GATK4_VQSR.out.versions)
    } else {
        GATK4_HARDFILTER(GATK4_MERGEVCFS_COHORT.out.vcf, ch_fasta, ch_fai, ch_dict)
        ch_filtered = GATK4_HARDFILTER.out.vcf
        ch_versions = ch_versions.mix(GATK4_HARDFILTER.out.versions)
    }

    BCFTOOLS_NORM_GERMLINE(ch_filtered, ch_fasta, ch_fai)
    BCFTOOLS_VIEW_PASS_GERMLINE(BCFTOOLS_NORM_GERMLINE.out.vcf)
    BCFTOOLS_STATS_GERMLINE(BCFTOOLS_NORM_GERMLINE.out.vcf, ch_fasta)
    ch_versions = ch_versions
        .mix(BCFTOOLS_NORM_GERMLINE.out.versions)
        .mix(BCFTOOLS_VIEW_PASS_GERMLINE.out.versions)
        .mix(BCFTOOLS_STATS_GERMLINE.out.versions)

    emit:
    vcf      = BCFTOOLS_NORM_GERMLINE.out.vcf       // [meta, vcf, tbi] all records, FILTER set
    pass_vcf = BCFTOOLS_VIEW_PASS_GERMLINE.out.vcf  // [meta, vcf, tbi] PASS only
    stats    = BCFTOOLS_STATS_GERMLINE.out.stats
    versions = ch_versions
}
