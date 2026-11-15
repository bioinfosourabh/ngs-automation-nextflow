# Output

Paths are relative to `--outdir`.

## Preprocessing

| Path | Content |
|---|---|
| `preprocessing/fastp/<sample>/<sample>_<lane>.fastp.{html,json}` | fastp reports per lane (trimmed FASTQ only with `--save_trimmed`) |
| `preprocessing/aligned/<sample>/` | lane-level sorted BAMs (`--save_aligned`) |
| `preprocessing/markduplicates/<sample>/<sample>.md.metrics.txt` | Picard duplicate metrics (BAM with `--save_markdup`) |
| `preprocessing/bqsr/<sample>/<sample>.recal.table` | BaseRecalibrator table |
| `preprocessing/recalibrated/<sample>/<sample>.recal.{bam,bai}` | analysis-ready BAM |
| `reference/` | FASTA index, dictionary and aligner index if built (`--save_reference`) |

## QC

| Path | Content |
|---|---|
| `qc/samtools/<sample>/` | `samtools stats` and `flagstat` |
| `qc/mosdepth/<sample>/` | coverage summary and distributions; per-target (`--intervals`) or 1 kb window depth |
| `qc/bcftools_stats/` | `bcftools stats` of the final call sets |
| `multiqc/multiqc_report.html` | combined report |

## Germline

`variant_calling/germline/joint_germline/`

| File | Content |
|---|---|
| `joint_germline.raw.vcf.gz` | joint-genotyped cohort VCF before filtering |
| `joint_germline.{hardfiltered,vqsr}.norm.vcf.gz` | filtered (FILTER column set), multiallelics split, indels left-aligned |
| `joint_germline.{hardfiltered,vqsr}.norm.pass.vcf.gz` | PASS records only |
| `variant_calling/germline/vqsr/` | recalibration and tranches files (`--germline_filter vqsr`) |
| `variant_calling/germline/<sample>/<sample>.g.vcf.gz` | per-sample GVCF (`--save_gvcf`) |

## Somatic

`variant_calling/somatic/<tumour>_vs_<normal>/` or `<tumour>_tumour_only/`

| File | Content |
|---|---|
| `*.mutect2.unfiltered.vcf.gz` | gathered Mutect2 calls |
| `*.contamination.table`, `*.segmentation.table` | CalculateContamination output (`--contamination_variants`) |
| `*.mutect2.filtered.vcf.gz` | FilterMutectCalls output |
| `*.mutect2.filteringStats.tsv` | filtering statistics |
| `*.mutect2.filtered.norm.vcf.gz` | normalised, all records |
| `*.mutect2.filtered.norm.pass.vcf.gz` | PASS records only |

## Annotation

`annotation/{germline,somatic}/<id>/*.vep.vcf.gz` and `*.vep_summary.html` (with `--run_vep`).

## Pipeline information

`pipeline_info/` holds the Nextflow execution report, timeline, trace, DAG and `software_versions.yml` (tool versions reported by each process).
