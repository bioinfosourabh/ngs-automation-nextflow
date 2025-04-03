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

## Contact
Sourabh Kumar ( Email: bioinfosourabh@gmail.com  )
