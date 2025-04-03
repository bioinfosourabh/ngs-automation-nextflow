#!/usr/bin/env nextflow

// Define input parameters
params.inputDir = "$baseDir"      // Assuming input directory is the directory from which the pipeline is launched
params.file_ref = "/media/prabudh/m1/hg38/hg38.fa"
params.known_site = "/media/prabudh/m1/vcf_file/dbsnp.vcf.gz"
params.platform = "ILLUMINA"
params.pon = "/media/prabudh/m1/vcf_file/somatic-hg38_1000g_pon.hg38.vcf"
params.af_only = "/media/prabudh/m1/vcf_file/somatic-hg38_af-only-gnomad.hg38.vcf"

// Output directories
outputDir = file("$params.inputDir/output")
dir_qc = file("$outputDir/fastqc")
dir_fastp = file("$outputDir/filtered_qc_report")
Dir_map = file("$outputDir/Mapsam")
Dir_vcf = file("$outputDir/Germline_VCF")

// Make sure output directories exist
process setupDirectories {
    output:
    file("${outputDir}/logs") into log_files

    script:
    """
    mkdir -p ${dir_qc}
    mkdir -p ${dir_fastp}
    mkdir -p ${Dir_map}
    mkdir -p ${Dir_vcf}
    """
}

// Define a channel for paired input files
Channel
    .fromFilePairs("${params.inputDir}/*_{1,2}.fastq.gz", size: 2)
    .set { sample_pairs }

// Step 1: Quality Control with Fastp
process fastp {
    input:
    tuple val(sample_name), file(fastq1), file(fastq2) from sample_pairs

    output:
    file("${dir_fastp}/${sample_name}_1_filtered.fastq.gz") into filtered_fastq

    maxForks 5  // Limit to 5 parallel executions

    script:
    """
    fastp -i ${fastq1} -I ${fastq2} --disable_length_filtering --qualified_quality_phred 30 \
          -o ${dir_fastp}/${sample_name}_1_filtered.fastq.gz \
          -O ${dir_fastp}/${sample_name}_2_filtered.fastq.gz \
          --html ${dir_fastp}/${sample_name}.html --thread 2
    """
}

// Step 2: Align with BWA
process bwa_mem {
    input:
    tuple val(sample_name), file(fastq1), file(fastq2) from filtered_fastq

    output:
    file("${Dir_map}/${sample_name}.bam") into bam_files

    maxForks 5  // Limit to 5 parallel executions

    script:
    """
    bwa mem -v 2 -M -t 32 -R "@RG\\tID:${sample_name}\\tPL:${params.platform}\\tLB:LIB_${sample_name}\\tSM:${sample_name}" \
        ${params.file_ref} ${fastq1} ${fastq2} | samtools sort -@ 12 -O BAM -o ${Dir_map}/${sample_name}.bam
    """
}

// Step 3: Mark Duplicates and Add Read Groups
process mark_duplicates {
    input:
    file(bam) from bam_files

    output:
    file("${bam.baseName}_markdup.bam") into dedup_bam

    maxForks 5  // Limit to 5 parallel executions

    script:
    """
    gatk MarkDuplicatesSpark -I ${bam} -O ${bam.baseName}_markdup.bam --spark-master local[12]
    """
}

// Step 4: Base Recalibration
process base_recalibrator {
    input:
    file(bam) from dedup_bam

    output:
    file("${bam.baseName}_recal.table") into recal_tables

    maxForks 5  // Limit to 5 parallel executions

    script:
    """
    gatk BaseRecalibrator --java-options "-Xmx4g -XX:ParallelGCThreads=2" -I ${bam} --known-sites ${params.known_site} \
        -O ${bam.baseName}_recal.table -R ${params.file_ref}
    """
}

// Step 5: Apply BQSR
process apply_bqsr {
    input:
    file(bam) from dedup_bam
    file(recal_table) from recal_tables

    output:
    file("${bam.baseName}_recal.bam") into recal_bam

    maxForks 5  // Limit to 5 parallel executions

    script:
    """
    gatk ApplyBQSR --java-options "-Xmx4g -XX:ParallelGCThreads=2" -bqsr ${recal_table} -I ${bam} \
        -O ${bam.baseName}_recal.bam -R ${params.file_ref}
    """
}

// Step 6: HaplotypeCaller
process haplotype_caller {
    input:
    file(bam) from recal_bam

    output:
    file("${bam.baseName}_detected_variants.vcf.gz") into gvcf_files

    maxForks 5  // Limit to 5 parallel executions

    script:
    """
    gatk HaplotypeCaller --java-options "-Xms20G -Xmx20G -XX:ParallelGCThreads=2" -R ${params.file_ref} \
        -I ${bam} -O ${bam.baseName}_detected_variants.vcf.gz -ERC GVCF
    """
}

// Step 7: GenotypeGVCFs
process genotype_gvcfs {
    input:
    file(gvcf) from gvcf_files

    output:
    file("${gvcf.baseName}_raw.vcf") into final_vcfs

    maxForks 5  // Limit to 5 parallel executions

    script:
    """
    gatk GenotypeGVCFs --java-options "-Xms2G -Xmx2g -XX:ParallelGCThreads=2" -R ${params.file_ref} \
        -V ${gvcf} -O ${gvcf.baseName}_raw.vcf --new-qual --max-alternate-alleles 2
    """
}

// Summary logs
workflow {
    setupDirectories()
    fastp()
    bwa_mem()
    mark_duplicates()
    base_recalibrator()
    apply_bqsr()
    haplotype_caller()
    genotype_gvcfs()
}


#### to run the pipeline 
nextflow run germline_variant_calling.nf -with-trace -with-report
