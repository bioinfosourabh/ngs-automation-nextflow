//
// Germline and somatic short-variant calling from paired-end FASTQ.
//
include { readSamplesheet               } from '../subworkflows/local/input_check'
include { PREPARE_GENOME                } from '../subworkflows/local/prepare_genome'
include { GERMLINE_JOINT                } from '../subworkflows/local/germline_joint'
include { SOMATIC_MUTECT2               } from '../subworkflows/local/somatic_mutect2'
include { FASTP                         } from '../modules/local/fastp'
include { BWAMEM2_MEM                   } from '../modules/local/bwamem2/mem'
include { BWA_MEM                       } from '../modules/local/bwa/mem'
include { GATK4_MARKDUPLICATES          } from '../modules/local/gatk4/markduplicates'
include { GATK4_BASERECALIBRATOR        } from '../modules/local/gatk4/baserecalibrator'
include { GATK4_APPLYBQSR               } from '../modules/local/gatk4/applybqsr'
include { SAMTOOLS_STATS                } from '../modules/local/samtools/stats'
include { MOSDEPTH                      } from '../modules/local/mosdepth'
include { ENSEMBLVEP_VEP                } from '../modules/local/ensemblvep/vep'
include { MULTIQC                       } from '../modules/local/multiqc'

// Optional single-file input with its tabix index; [] when not given.
def optionalVcf(path) {
    return path ? [file(path, checkIfExists: true), file("${path}.tbi", checkIfExists: true)] : [[], []]
}

workflow NGS_VARIANTS {
    main:
    ch_versions = channel.empty()
    ch_mqc      = channel.empty()

    def samples  = readSamplesheet(params.input, params.seq_platform)
    def statuses = samples.collect { s -> s[0].status }
    def has_normals = statuses.contains(0)
    def has_tumours = statuses.contains(1)
    ch_reads = channel.fromList(samples)

    //
    // Reference files
    //
    PREPARE_GENOME()
    ch_fasta = PREPARE_GENOME.out.fasta
    ch_fai   = PREPARE_GENOME.out.fai
    ch_dict  = PREPARE_GENOME.out.dict
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    def dbsnp = optionalVcf(params.dbsnp)
    def known_indels = params.known_indels ? params.known_indels.tokenize(',').collect { p -> optionalVcf(p.trim()) } : []
    ch_known_sites     = channel.value([dbsnp[0]] + known_indels.collect { k -> k[0] })
    ch_known_sites_tbi = channel.value([dbsnp[1]] + known_indels.collect { k -> k[1] })

    //
    // Read QC and trimming, alignment per lane
    //
    FASTP(ch_reads)
    ch_mqc      = ch_mqc.mix(FASTP.out.json.map { _meta, json -> json })
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    ch_lane_bam = channel.empty()
    if (params.aligner == 'bwa-mem2') {
        BWAMEM2_MEM(FASTP.out.reads, PREPARE_GENOME.out.index)
        ch_lane_bam = BWAMEM2_MEM.out.bam
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
    } else {
        BWA_MEM(FASTP.out.reads, PREPARE_GENOME.out.index)
        ch_lane_bam = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    }

    //
    // Merge lanes + mark duplicates, then BQSR
    //
    ch_md_in = ch_lane_bam
        .map { meta, bam ->
            def smeta = meta.findAll { k, _v -> !(k in ['id', 'lane', 'read_group', 'num_lanes']) } + [id: meta.sample]
            [groupKey(smeta, meta.num_lanes), bam]
        }
        .groupTuple()
        .map { key, bams -> [key.getGroupTarget(), bams.sort { b -> b.name }] }
    GATK4_MARKDUPLICATES(ch_md_in)
    ch_mqc      = ch_mqc.mix(GATK4_MARKDUPLICATES.out.metrics.map { _meta, m -> m })
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())

    GATK4_BASERECALIBRATOR(
        GATK4_MARKDUPLICATES.out.bam, ch_fasta, ch_fai, ch_dict,
        ch_known_sites, ch_known_sites_tbi, PREPARE_GENOME.out.intervals
    )
    ch_mqc      = ch_mqc.mix(GATK4_BASERECALIBRATOR.out.table.map { _meta, t -> t })
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first())

    GATK4_APPLYBQSR(
        GATK4_MARKDUPLICATES.out.bam.join(GATK4_BASERECALIBRATOR.out.table),
        ch_fasta, ch_fai, ch_dict
    )
    ch_bam = GATK4_APPLYBQSR.out.bam
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions.first())

    //
    // Alignment QC
    //
    SAMTOOLS_STATS(ch_bam, ch_fasta)
    MOSDEPTH(ch_bam, PREPARE_GENOME.out.intervals)
    ch_mqc = ch_mqc
        .mix(SAMTOOLS_STATS.out.stats.map { _meta, f -> f })
        .mix(SAMTOOLS_STATS.out.flagstat.map { _meta, f -> f })
        .mix(MOSDEPTH.out.global_txt.map { _meta, f -> f })
        .mix(MOSDEPTH.out.summary_txt.map { _meta, f -> f })
        .mix(MOSDEPTH.out.regions_txt.map { _meta, f -> f })
    ch_versions = ch_versions
        .mix(SAMTOOLS_STATS.out.versions.first())
        .mix(MOSDEPTH.out.versions.first())

    ch_to_annotate = channel.empty()

    //
    // Germline joint calling on normal samples
    //
    if (params.run_germline && !has_normals) {
        log.warn("--run_germline is set but the samplesheet has no normal (status 0) samples; skipping germline calling")
    }
    if (params.run_germline && has_normals) {
        def vqsr_names = [:]
        def vqsr_files = []
        if (params.germline_filter == 'vqsr') {
            [hapmap: params.vqsr_hapmap, omni: params.vqsr_omni, g1k: params.vqsr_1000g,
             mills: params.vqsr_mills, dbsnp: params.dbsnp, axiom: params.vqsr_axiom].each { name, path ->
                if (path) {
                    def pair = optionalVcf(path)
                    vqsr_files.addAll(pair)
                    vqsr_names[name] = pair[0].name
                }
            }
        }
        GERMLINE_JOINT(
            ch_bam.filter { meta, _bam, _bai -> meta.status == 0 },
            PREPARE_GENOME.out.scatter,
            ch_fasta, ch_fai, ch_dict,
            dbsnp[0], dbsnp[1],
            channel.value(vqsr_files),
            vqsr_names
        )
        ch_mqc         = ch_mqc.mix(GERMLINE_JOINT.out.stats.map { _meta, f -> f })
        ch_to_annotate = ch_to_annotate.mix(GERMLINE_JOINT.out.pass_vcf)
        ch_versions    = ch_versions.mix(GERMLINE_JOINT.out.versions)
    }

    //
    // Somatic calling for tumour samples
    //
    if (params.run_somatic && !has_tumours) {
        log.warn("--run_somatic is set but the samplesheet has no tumour (status 1) samples; skipping somatic calling")
    }
    if (params.run_somatic && has_tumours) {
        def germline_resource = optionalVcf(params.germline_resource)
        def pon    = optionalVcf(params.pon)
        def contam = optionalVcf(params.contamination_variants)
        SOMATIC_MUTECT2(
            ch_bam,
            PREPARE_GENOME.out.scatter,
            ch_fasta, ch_fai, ch_dict,
            PREPARE_GENOME.out.intervals,
            germline_resource[0], germline_resource[1],
            pon[0], pon[1],
            contam[0], contam[1]
        )
        ch_mqc         = ch_mqc.mix(SOMATIC_MUTECT2.out.stats.map { _meta, f -> f })
        ch_to_annotate = ch_to_annotate.mix(SOMATIC_MUTECT2.out.pass_vcf)
        ch_versions    = ch_versions.mix(SOMATIC_MUTECT2.out.versions)
    }

    //
    // Optional annotation of PASS variants with Ensembl VEP
    //
    if (params.run_vep) {
        ENSEMBLVEP_VEP(ch_to_annotate, ch_fasta, channel.value(file(params.vep_cache, checkIfExists: true)))
        ch_mqc      = ch_mqc.mix(ENSEMBLVEP_VEP.out.report.map { _meta, f -> f })
        ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())
    }

    //
    // Software versions and MultiQC
    //
    ch_versions_yml = ch_versions
        .unique { f -> f.text }
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true, newLine: true)

    ch_report = channel.empty()
    if (!params.skip_multiqc) {
        ch_multiqc_config = channel.value(file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true))
        MULTIQC(ch_mqc.collect(), ch_multiqc_config)
        ch_report = MULTIQC.out.report
    }

    emit:
    bam            = ch_bam
    versions       = ch_versions_yml
    multiqc_report = ch_report
}
