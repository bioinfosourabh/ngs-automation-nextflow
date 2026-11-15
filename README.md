# NGS Analysis Automation with Nextflow: Germline Variant Calling

A prototype Nextflow implementation of a germline short-variant calling workflow (fastp, BWA-MEM, GATK4) for paired-end WES/WGS data on hg38. It ports the per-sample steps of the Bash [Germline-Variant-Calling-Pipeline](https://github.com/bioinfosourabh/Germline-Variant-Calling-Pipeline) into Nextflow processes so that samples can run in parallel.

**Status: prototype.** The script uses DSL1-style channel syntax (`from` / `into`) together with a DSL2 `workflow {}` block. It needs to be ported fully to DSL2 before it will run on current Nextflow releases.

## Workflow

Script: [`scripts/germline_variant_calling.nf`](scripts/germline_variant_calling.nf)

| Process | Tool | Description |
|---|---|---|
| `setupDirectories` | bash | creates the output folders |
| `fastp` | fastp | read trimming and filtering (Phred 30 qualified-base threshold), HTML report |
| `bwa_mem` | BWA-MEM, samtools | alignment with read groups, coordinate-sorted BAM |
| `mark_duplicates` | GATK MarkDuplicatesSpark | duplicate marking |
| `base_recalibrator` | GATK BaseRecalibrator | BQSR table from dbSNP known sites |
| `apply_bqsr` | GATK ApplyBQSR | recalibrated BAM |
| `haplotype_caller` | GATK HaplotypeCaller | per-sample GVCF (`-ERC GVCF`) |
| `genotype_gvcfs` | GATK GenotypeGVCFs | genotypes each sample's GVCF into a VCF |

Each process is limited to 5 parallel tasks (`maxForks 5`). Genotyping is done per sample; there is no CombineGVCFs or joint-genotyping step yet.

## Input

- Paired-end FASTQ files named `<sample>_1.fastq.gz` and `<sample>_2.fastq.gz` in `params.inputDir` (default: the directory containing the script)
- hg38 reference FASTA with BWA and samtools indexes and a sequence dictionary
- dbSNP VCF for BQSR

## Parameters

Set these in the `params` block at the top of the script (or on the command line, e.g. `--file_ref /path/to/hg38.fa`):

```groovy
params.inputDir   = "$baseDir"
params.file_ref   = "/path/to/hg38.fa"
params.known_site = "/path/to/dbsnp.vcf.gz"
params.platform   = "ILLUMINA"
```

`params.pon` and `params.af_only` (panel of normals and gnomAD allele-frequency VCF) are also defined but are not used by any process; they are left over from the somatic pipeline.

## Usage

```sh
nextflow run scripts/germline_variant_calling.nf -with-trace -with-report -resume
```

Useful extra options: `-with-dag dag.png` (workflow graph) and `-with-timeline timeline.html` (run timeline). No Docker or Conda profile is defined yet, so the tools must be on `PATH`.

## Output

```
output/
├── filtered_qc_report/   fastp-filtered FASTQ files and HTML reports
├── Mapsam/               sorted BAM files
├── fastqc/               created but not used
└── Germline_VCF/         created but not used
```

Duplicate-marked and recalibrated BAMs, GVCFs and the final `*_raw.vcf` files are written to the Nextflow `work/` directories, because no `publishDir` is set.

## Known gaps

- Not yet ported to Nextflow DSL2 (see Status above)
- No joint genotyping across samples
- No `publishDir`, `nextflow.config` or container profiles

## License

MIT, see [LICENSE](LICENSE).

## Author

Sourabh Kumar · [bioinfosourabh@gmail.com](mailto:bioinfosourabh@gmail.com) · [bioinfosourabh.netlify.app](https://bioinfosourabh.netlify.app) · ORCID [0009-0009-5840-9830](https://orcid.org/0009-0009-5840-9830)
