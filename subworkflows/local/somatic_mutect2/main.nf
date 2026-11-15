//
// Somatic short variants with Mutect2 (tumour-normal or tumour-only):
// scattered calling, gather, read-orientation model, contamination
// estimation, FilterMutectCalls and normalisation.
//
include { GATK4_MUTECT2                                    } from '../../../modules/local/gatk4/mutect2'
include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_MUTECT2       } from '../../../modules/local/gatk4/mergevcfs'
include { GATK4_MERGEMUTECTSTATS                           } from '../../../modules/local/gatk4/mergemutectstats'
include { GATK4_LEARNREADORIENTATIONMODEL                  } from '../../../modules/local/gatk4/learnreadorientationmodel'
include { GATK4_GETPILEUPSUMMARIES                         } from '../../../modules/local/gatk4/getpileupsummaries'
include { GATK4_CALCULATECONTAMINATION                     } from '../../../modules/local/gatk4/calculatecontamination'
include { GATK4_FILTERMUTECTCALLS                          } from '../../../modules/local/gatk4/filtermutectcalls'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_SOMATIC           } from '../../../modules/local/bcftools/norm'
include { BCFTOOLS_VIEW_PASS as BCFTOOLS_VIEW_PASS_SOMATIC } from '../../../modules/local/bcftools/view'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_SOMATIC         } from '../../../modules/local/bcftools/stats'

workflow SOMATIC_MUTECT2 {
    take:
    ch_bam                  // [meta, bam, bai] all samples
    ch_scatter              // [chunk, interval_list]
    ch_fasta
    ch_fai
    ch_dict
    ch_intervals            // BED or []
    ch_germline_resource
    ch_germline_resource_tbi
    ch_pon                  // VCF or []
    ch_pon_tbi
    ch_contam_vcf           // VCF or [] (contamination estimation is skipped when empty)
    ch_contam_tbi

    main:
    ch_versions = channel.empty()

    ch_tumour = ch_bam.filter { meta, _bam, _bai -> meta.status == 1 }
    ch_normal = ch_bam
        .filter { meta, _bam, _bai -> meta.status == 0 }
        .map { meta, bam, bai -> [meta.sample, meta, bam, bai] }

    // Tumour-normal pairs (matched by patient) and tumour-only samples.
    ch_paired = ch_tumour
        .filter { meta, _bam, _bai -> meta.normal_id }
        .map { meta, bam, bai -> [meta.normal_id, meta, bam, bai] }
        .combine(ch_normal, by: 0)
        .map { _normal_id, tmeta, tbam, tbai, nmeta, nbam, nbai ->
            def pair = [id: "${tmeta.sample}_vs_${nmeta.sample}", patient: tmeta.patient, tumour_id: tmeta.sample, normal_id: nmeta.sample, analysis: 'somatic']
            [pair, [tbam, nbam], [tbai, nbai]]
        }
    ch_tumour_only = ch_tumour
        .filter { meta, _bam, _bai -> !meta.normal_id }
        .map { meta, bam, bai ->
            def single = [id: "${meta.sample}_tumour_only", patient: meta.patient, tumour_id: meta.sample, normal_id: null, analysis: 'somatic']
            [single, [bam], [bai]]
        }
    ch_pairs = ch_paired.mix(ch_tumour_only)

    // Scattered calling
    ch_m2_in = ch_pairs
        .combine(ch_scatter)
        .map { meta, bams, bais, chunk, interval -> [meta + [chunk: chunk], bams, bais, interval] }
    GATK4_MUTECT2(ch_m2_in, ch_fasta, ch_fai, ch_dict, ch_germline_resource, ch_germline_resource_tbi, ch_pon, ch_pon_tbi)
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

    // Gather per pair
    ch_vcf_parts = GATK4_MUTECT2.out.vcf
        .map { meta, vcf, _tbi -> [meta.findAll { k, _v -> k != 'chunk' }, vcf] }
        .groupTuple()
        .map { meta, vcfs -> [meta, vcfs.sort { v -> v.name }] }
    GATK4_MERGEVCFS_MUTECT2(ch_vcf_parts, ch_dict, '.mutect2.unfiltered')

    ch_stats_parts = GATK4_MUTECT2.out.stats
        .map { meta, stats -> [meta.findAll { k, _v -> k != 'chunk' }, stats] }
        .groupTuple()
        .map { meta, stats -> [meta, stats.sort { s -> s.name }] }
    GATK4_MERGEMUTECTSTATS(ch_stats_parts)

    ch_f1r2_parts = GATK4_MUTECT2.out.f1r2
        .map { meta, f1r2 -> [meta.findAll { k, _v -> k != 'chunk' }, f1r2] }
        .groupTuple()
        .map { meta, f1r2 -> [meta, f1r2.sort { f -> f.name }] }
    GATK4_LEARNREADORIENTATIONMODEL(ch_f1r2_parts)

    ch_versions = ch_versions
        .mix(GATK4_MERGEVCFS_MUTECT2.out.versions.first())
        .mix(GATK4_MERGEMUTECTSTATS.out.versions.first())
        .mix(GATK4_LEARNREADORIENTATIONMODEL.out.versions.first())

    ch_filter_in = GATK4_MERGEVCFS_MUTECT2.out.vcf
        .join(GATK4_MERGEMUTECTSTATS.out.stats)
        .join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior)

    if (params.contamination_variants) {
        // Pileups at common germline sites for every tumour and every matched normal.
        ch_pileup_in = ch_bam.filter { meta, _bam, _bai -> meta.status == 1 || meta.has_tumour }
        GATK4_GETPILEUPSUMMARIES(ch_pileup_in, ch_fasta, ch_fai, ch_dict, ch_contam_vcf, ch_contam_tbi, ch_intervals)
        ch_tables = GATK4_GETPILEUPSUMMARIES.out.table.map { meta, table -> [meta.sample, table] }

        ch_cc_base = ch_pairs
            .map { meta, _bams, _bais -> [meta.tumour_id, meta] }
            .combine(ch_tables, by: 0)
            .map { _tumour_id, meta, ttable -> [meta.normal_id ?: 'NO_NORMAL', meta, ttable] }
        ch_cc_paired = ch_cc_base
            .filter { key, _meta, _ttable -> key != 'NO_NORMAL' }
            .combine(ch_tables, by: 0)
            .map { _key, meta, ttable, ntable -> [meta, ttable, ntable] }
        ch_cc_tonly = ch_cc_base
            .filter { key, _meta, _ttable -> key == 'NO_NORMAL' }
            .map { _key, meta, ttable -> [meta, ttable, []] }

        GATK4_CALCULATECONTAMINATION(ch_cc_paired.mix(ch_cc_tonly))
        ch_filter_in = ch_filter_in.join(GATK4_CALCULATECONTAMINATION.out.tables)
        ch_versions = ch_versions
            .mix(GATK4_GETPILEUPSUMMARIES.out.versions.first())
            .mix(GATK4_CALCULATECONTAMINATION.out.versions.first())
    } else {
        ch_filter_in = ch_filter_in.map { meta, vcf, tbi, stats, prior -> [meta, vcf, tbi, stats, prior, [], []] }
    }

    GATK4_FILTERMUTECTCALLS(ch_filter_in, ch_fasta, ch_fai, ch_dict)
    BCFTOOLS_NORM_SOMATIC(GATK4_FILTERMUTECTCALLS.out.vcf, ch_fasta, ch_fai)
    BCFTOOLS_VIEW_PASS_SOMATIC(BCFTOOLS_NORM_SOMATIC.out.vcf)
    BCFTOOLS_STATS_SOMATIC(BCFTOOLS_VIEW_PASS_SOMATIC.out.vcf, ch_fasta)
    ch_versions = ch_versions
        .mix(GATK4_FILTERMUTECTCALLS.out.versions.first())
        .mix(BCFTOOLS_NORM_SOMATIC.out.versions.first())
        .mix(BCFTOOLS_VIEW_PASS_SOMATIC.out.versions.first())
        .mix(BCFTOOLS_STATS_SOMATIC.out.versions.first())

    emit:
    vcf      = BCFTOOLS_NORM_SOMATIC.out.vcf        // [meta, vcf, tbi] all records, FILTER set
    pass_vcf = BCFTOOLS_VIEW_PASS_SOMATIC.out.vcf   // [meta, vcf, tbi] PASS only
    stats    = BCFTOOLS_STATS_SOMATIC.out.stats
    versions = ch_versions
}
