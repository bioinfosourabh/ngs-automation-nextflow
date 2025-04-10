# NGS Analysis Automation using NextFlow
## Germline Variant Calling Workflow (Nextflow)
Overview

This repository implements a modular and reproducible germline variant calling pipeline using Nextflow. It automates a standard short-read variant discovery workflow based on the GATK Best Practices, and is designed for analyzing large-scale whole exome or whole genome sequencing data.
The pipeline supports parallel execution, reproducibility, and easy scaling across local machines, HPC clusters, or cloud environments.

## Features
- Developed with Nextflow DSL
- Modular and customizable
- Compatible with GATK 4.x, BWA, Samtools, fastp
- Supports per-sample GVCF generation and joint genotyping
- Produces structured outputs and summary logs

## Input Requirements
To run the pipeline, you will need:
- Paired-end FASTQ files (*_1.fastq.gz and *_2.fastq.gz)
- Reference genome FASTA file (hg38.fa) and index files
- dbSNP VCF file (known sites)
- gnomAD allele frequency VCF
- Panel of Normals (PoN) for filtering technical artifacts

## Output Structure
The pipeline generates the following output directories:
```
output/
├── fastqc/                  # Raw QC reports
├── filtered_qc_report/      # Trimmed/filtered FASTQs (via fastp)
├── Mapsam/                  # Aligned and sorted BAM files
├── Germline_VCF/            # GVCF and raw VCF files
└── logs/                    # Log files for each Nextflow process
```

## Pipeline Steps
| Step  | Tool  | Description   |
|----------|----------|----------|
| Directory setup	  | bash  | 	Initializes output folders   |
| Quality filtering	  | fastp  | 	Adapter trimming and QC   |
| Alignment	  | BWA-MEM  | 	Aligns reads to reference   |
| Duplicate marking  | GATK MarkDuplicatesSpark  | 	Marks PCR duplicates   |
| Base recalibration  | GATK BaseRecalibrator  | 	Uses known SNPs for BQSR   |
| Apply recalibration	  | GATK ApplyBQSR  | 	Applies BQSR table to BAM   |
| Variant calling  | GATK HaplotypeCaller  | 	GVCF mode per sample   |
| Joint genotyping  | GATK GenotypeGVCFs  | 	Combines GVCFs into VCF   |

## Nextflow Workflow
The Nextflow script defines the following structure:
```
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
```
Each process block (e.g., fastp, bwa_mem, mark_duplicates, etc.) is defined in the .nf script and configured with inputs, outputs, and commands for execution. Channels are used to handle paired-end files and pass them between processes.

## Running the Pipeline
Make sure you have Nextflow installed and properly configured.
Run the pipeline with:
```
nextflow run germline_variant_calling.nf -with-trace -with-report -resume
```
Add optional flags as needed:
-with-dag dag.png to visualize the pipeline structure
-with-timeline timeline.html to assess performance
-profile docker for containerized execution (if implemented)


## Parameter Configuration
Inside the germline_variant_calling.nf script, edit the params block to match your file paths:
```
params.inputDir     = "$baseDir"
params.file_ref     = "/path/to/hg38.fa"
params.known_site   = "/path/to/dbsnp.vcf.gz"
params.pon          = "/path/to/pon.vcf.gz"
params.af_only      = "/path/to/gnomad.af.vcf.gz"
params.platform     = "ILLUMINA"
```

## Whole Nextflow Pipeline
```
#!/usr/bin/env nextflow

// Define input parameters
params.inputDir = "$baseDir"      // Assuming input directory is the directory from which the pipeline is launched
params.file_ref = "/path/to/hg38.fa"
params.known_site = "/path/to/dbsnp.vcf.gz"
params.platform = "ILLUMINA"
params.pon = "/path/to/somatic-hg38_1000g_pon.hg38.vcf"
params.af_only = "/path/to/somatic-hg38_af-only-gnomad.hg38.vcf"

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
```

## Contact
Sourabh Kumar ( Email: bioinfosourabh@gmail.com  )
