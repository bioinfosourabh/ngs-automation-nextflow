# Usage

## Samplesheet

| Column | Required | Description |
|---|---|---|
| `patient` | yes | Patient or individual ID. Tumours are paired with the normal of the same patient. |
| `sample` | yes | Sample ID, used as the read-group `SM` tag and in output names. Letters, digits, `.`, `_`, `-`. |
| `status` | no | `0` normal (default), `1` tumour. |
| `lane` | no | Lane or library ID. If missing, rows of the same sample are numbered 1, 2, ... |
| `fastq_1`, `fastq_2` | yes | Gzipped paired FASTQ files (`.fastq.gz` or `.fq.gz`). Absolute paths, URLs/S3 URIs, or paths relative to the samplesheet's directory. |

Checks made before the run starts: required columns present; sample and patient IDs well formed; status is 0 or 1; FASTQ files exist; each sample/lane pair listed once; patient and status consistent across a sample's rows; at most one normal per patient.

## Analysis modes

| Samples in the sheet | What runs |
|---|---|
| normals only | germline joint genotyping |
| tumours with a normal of the same patient | Mutect2 tumour-normal, and germline calling of the normals |
| tumours without a normal | Mutect2 tumour-only |

Turn either branch off with `--run_germline false` or `--run_somatic false`.

## Reference files

For GRCh38, the GATK resource bundle provides everything needed:

| Parameter | File |
|---|---|
| `--fasta` | `Homo_sapiens_assembly38.fasta` (with `.fai` and `.dict`, or let the pipeline build them) |
| `--bwa_index` | directory with a pre-built index; building a BWA-MEM2 index for GRCh38 needs about 64 GB of memory |
| `--dbsnp` | `Homo_sapiens_assembly38.dbsnp138.vcf.gz` |
| `--known_indels` | `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz,Homo_sapiens_assembly38.known_indels.vcf.gz` |
| `--germline_resource` | `af-only-gnomad.hg38.vcf.gz` |
| `--pon` | `1000g_pon.hg38.vcf.gz`, or a panel of normals built from your own normals with the same assay |
| `--contamination_variants` | `small_exac_common_3.hg38.vcf.gz` |
| `--vqsr_hapmap`, `--vqsr_omni`, `--vqsr_1000g`, `--vqsr_mills`, `--vqsr_axiom` | `hapmap_3.3.hg38.vcf.gz`, `1000G_omni2.5.hg38.vcf.gz`, `1000G_phase1.snps.high_confidence.hg38.vcf.gz`, `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`, `Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz` |

All VCFs must be bgzipped and have a `.tbi` index next to them.

## Exome and panel data

Pass the capture targets as a BED file with `--intervals`. They restrict BQSR, variant calling and contamination estimation, and mosdepth reports coverage over them. GATK recommends `--interval_padding 100` for exomes. For VQSR on exomes, remove the `DP` annotation (it is not used by default here) and consider `--germline_filter hard` for small cohorts.

## Scatter-gather

The calling regions are split into `--scatter_count` chunks with GATK SplitIntervals (whole contigs are not split). HaplotypeCaller, GenomicsDBImport/CombineGVCFs, GenotypeGVCFs and Mutect2 run once per chunk, and the results are merged. Raise the count for whole genomes on a cluster or cloud (for example 24 to 50).

## Customising tools

Extra arguments can be passed to any process with `ext.args` in a config file given with `-c`:

```groovy
process {
    withName: 'GATK4_MUTECT2' {
        ext.args = '--genotype-germline-sites true --genotype-pon-sites true'
    }
    withName: 'GATK4_VQSR' {
        ext.snp_args   = '-an DP'
        ext.indel_args = '-an DP'
    }
}
```

## Resources

Defaults are set per process label in `conf/base.config` and capped by `process.resourceLimits` (16 CPUs, 64 GB, 48 h). For a workstation with fewer resources:

```groovy
process.resourceLimits = [cpus: 8, memory: 32.GB, time: 24.h]
```

## Resuming

Add `-resume` to reuse completed tasks after a failure or a parameter change that affects only later steps. Inputs to multi-file steps are sorted by name so that cached results are found again.

## Running on AWS Batch

```bash
nextflow run bioinfosourabh/ngs-automation-nextflow -profile awsbatch \
    --awsqueue my-queue --awsregion eu-west-2 \
    -work-dir s3://my-bucket/work --outdir s3://my-bucket/results \
    --input s3://my-bucket/samplesheet.csv ...
```

Infrastructure templates and a worked example are in [aws-genomics-workflow](https://github.com/bioinfosourabh/aws-genomics-workflow).
