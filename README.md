# Whole Exome Sequencing Analysis Pipeline

A bioinformatics pipeline for processing whole exome sequencing data aligned to the hg38/GRCh38 reference genome.

## Pipeline Overview

This pipeline performs the following steps:
1. Quality control of raw FASTQ files (FastQC)
2. Adapter and Quality trimming and quality filtering (BBTools)
3. Alignment to hg38 reference genome (BWA-MEM)
4. SAM to BAM conversion, sorting, and indexing (samtools)
5. Duplicate marking and removal (Picard)
6. Local realignment around indels (GATK)
7. Base quality score recalibration (GATK)
8. Variant calling (GATK HaplotypeCaller)
9. Variant filtering (GATK VariantFiltration)

## Requirements

### Software Dependencies
- BWA (v0.7.17 or later)
- samtools (v1.10 or later)
- Picard Tools (v2.23 or later)
- GATK (v3.8 or later)
- BBTools (for bbduk.sh)
- FastQC (v0.11.9 or later)
- Java (v1.8)

### Reference Files Needed
- hg38/GRCh38 reference genome (BWA indexed)
- Known variant sites:
  - Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
  - 1000G_phase1.snps.high_confidence.hg38.vcf.gz
  - dbsnp154_hg38.vcf.gz
- Target capture regions BED file (optional)

## Installation

Clone this repository:

git clone https://github.com/Maroua2906/Wes_Analysis_Pipeline.git

### Usage
bash scripts/wes_pipeline.sh

### Customization
**Before running:**
1. Replace all `/PATH/` placeholders in the script with your actual paths
2. Update sample-specific information:
   - `@RG` header in BWA command (ID, SM, DT, LB, PL, CN)
3. Adjust resource parameters:
   - Threads (`-t` in BWA, `-nct` in GATK)
   - Memory (`-Xmx12g` in Java tools)

### Input Requirements
- **Paired-end FASTQ files** (gzipped)
- **Naming convention**:  
  `{SAMPLE}_1.fastq.gz` and `{SAMPLE}_2.fastq.gz`

## Output Files
Directory structure created:
```
results/
├── quality_control/       # FastQC reports + trimmed FASTQs
├── alignment/             # Intermediate SAM files
├── make_bam/              # Initial BAMs
├── process_bam/           # Processed BAMs (deduped, realigned)
└── variants_calls/        # Final VCFs
```

## Important Notes
- **Resource recommendations**:
  - Minimum 12GB RAM (adjust `-Xmx` in Java commands)
  - Multi-threading: Set `-t` (BWA) and `-nct` (GATK) 
- **QC aggregation**:  
  Consider adding MultiQC (`multiqc .`) in quality_control/ directories

## Troubleshooting

- Path errors : Double-check all reference file paths
- Java memory errors : Increase `-Xmx` (e.g., `-Xmx24g`)
- Threading crashes : Reduce `-t`/`-nct` values
- Empty output : Verify FASTQ inputs exist and are not corrupted 

