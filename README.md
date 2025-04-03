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

Input Requirements
To run the pipeline, you will need:
- Paired-end FASTQ files (*_1.fastq.gz and *_2.fastq.gz)
- Reference genome FASTA file (hg38.fa) and index files
- dbSNP VCF file (known sites)
- gnomAD allele frequency VCF
- Panel of Normals (PoN) for filtering technical artifacts

Output Structure
The pipeline generates the following output directories:
output/
├── fastqc/                  # Raw QC reports
├── filtered_qc_report/      # Trimmed/filtered FASTQs (via fastp)
├── Mapsam/                  # Aligned and sorted BAM files
├── Germline_VCF/            # GVCF and raw VCF files
└── logs/                    # Log files for each Nextflow process



  
