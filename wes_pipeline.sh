#!/bin/bash
###############################################################################
# Whole Exome Sequencing Analysis Pipeline
# 
# Description: 
#   This script performs end-to-end analysis of exome sequencing data aligned to
#   hg38/GRCh38, including QC, alignment, variant calling, and filtering.
#
# Key Features:
#   - Processes paired-end FASTQ files
#   - Performs adapter trimming and quality control
#   - Aligns to hg38 using BWA-MEM
#   - Uses GATK best practices for variant calling


# Usage:
#   bash wes_pipeline.sh
#
# Requirements: BWA, samtools, Picard, GATK, BBTools, FastQC
#

### INITIALIZATION #############################################################

echo "Starting WES analysis pipeline - $(date)"

## Create directory structure
echo "Creating output directories..."
mkdir -p /PATH/exome_analysis_results_hg38/Test1_wes_analysis/{Test1_quality_control,Test1_alignment,Test1_make_bam,Test1_process_bam,Test1_variants_calls}

### QUALITY CONTROL ############################################################

echo "=== QUALITY CONTROL PHASE ==="

## Raw FASTQ QC
echo "Running FastQC on raw reads..."
fastqc /PATH/Test1_1.fastq.gz --outdir /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control
fastqc /PATH/Test1_2.fastq.gz --outdir /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control

## Adapter trimming
echo "Trimming adapters with bbduk..."
/usr/local/bbmap/bbduk.sh \
  in1=/PATH/Test1_1.fastq.gz \
  in2=/PATH/Test1_2.fastq.gz \
  out=/PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control/Test1_adapters_trimmed#.fastq.gz \
  ref=/usr/local/bbmap/resources/adapters.fa \
  ktrim=r k=23 mink=11 t=12

## Trimmed FASTQ QC
echo "Running FastQC on trimmed reads..."
fastqc /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control/Test1_adapters_trimmed1.fastq.gz --outdir /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control
fastqc /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control/Test1_adapters_trimmed2.fastq.gz --outdir /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control

## Remove short reads
echo "Filtering short reads (<30bp)..."
/usr/local/bbmap/bbduk.sh \
  in1=/PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control/Test1_adapters_trimmed1.fastq.gz \
  in2=/PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control/Test1_adapters_trimmed2.fastq.gz \
  out=/PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control/Test1_adapters_trimmed_l30_#.fastq.gz \
  minlen=30

### ALIGNMENT ##################################################################

echo "=== ALIGNMENT PHASE ==="

echo "Aligning to hg38 with BWA-MEM..."
bwa mem -M -t 10 \
  -R '@RG\tID:Test1\tSM:Test1\tDT:DATE\tLB:XXXXXXX\tPL:ILLUMINA\tCN:XXXXXXX' \
  /PATH/BWAIndex/genome.fa \
  /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control/Test1_adapters_trimmed_l30_1.fastq.gz \
  /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_quality_control/Test1_adapters_trimmed_l30_2.fastq.gz \
  > /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_alignment/Test1.sam

### BAM PROCESSING #############################################################

echo "=== BAM PROCESSING PHASE ==="

## SAM to BAM conversion
echo "Converting SAM to BAM..."
samtools view -bt /PATH/genome.fa \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_make_bam/Test1.bam \
  /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_alignment/Test1.sam

## Sort and index
echo "Sorting and indexing BAM..."
samtools sort /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_make_bam/Test1.bam \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_make_bam/Test1.srt.bam
samtools index /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_make_bam/Test1.srt.bam

## Picard processing
echo "Running Picard tools..."

## QC of the BAM file
java -Xmx12g -jar /usr/local/picard.jar ValidateSamFile \
  I=/PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_make_bam/Test1.srt.bam \
  MODE=SUMMARY

java -Xmx12g -jar /usr/local/picard.jar FixMateInformation \
  I=/PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_make_bam/Test1.srt.bam \
  O=/PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.bam \
  VALIDATION_STRINGENCY=LENIENT

samtools sort /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.bam \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.bam
samtools index /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.bam

## Mark duplicates
  echo "Marking duplicates..."
java -Xmx12g -jar /usr/local/picard.jar MarkDuplicates \
  I=/PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.bam \
  O=/PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.bam \
  M=/PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1_marked_dup_metrics.txt \
  VALIDATION_STRINGENCY=LENIENT \
  CREATE_INDEX=true

### QC Post Alignement ############################################################

echo "=== VARIANT CALLING PHASE ==="

## Local realignment
echo "Performing local realignment..."
#Identify what regions need to be realigned
/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Xmx12g -jar /usr/local/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R /PATH/genome.fa \
  -I /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.bam \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.intervals \
  -L /PATH/Exome_RefSeq_targets_hg38.bed \
  --known /PATH/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known /PATH/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  -nt 10
#Perform the local realignment 
/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Xmx12g -jar /usr/local/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R /PATH/genome.fa \
  -I /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.bam \
  -targetIntervals /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.intervals \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.realign.bam \
  -L /PATH/Exome_RefSeq_targets_hg38.bed \
  -known /PATH/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -known /PATH/1000G_phase1.snps.high_confidence.hg38.vcf.gz

## Base quality recalibration
echo "Running base quality recalibration..."
/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Xmx12g -jar /usr/local/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -R /PATH/genome.fa \
  -I /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.realign.bam \
  -L /PATH/Exome_RefSeq_targets_hg38.bed \
  --knownSites /scratch/dbsnp_files_hg38/dbsnp154_hg38.vcf.gz \
  --knownSites /PATH/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --knownSites /PATH/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  -nct 12 \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.recal_data.table \
  -l INFO

/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Xmx12g -jar /usr/local/GenomeAnalysisTK.jar \
  -T PrintReads \
  -R /PATH/genome.fa \
  -I /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.realign.bam \
  -BQSR /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.recal_data.table \
  -nct 12 \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.realign.recal.bam \
  -l INFO

## Final BAM processing
echo "Creating final BAM..."
samtools sort /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.realign.recal.bam \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.realign.recal.srt.bam
samtools index /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.realign.recal.srt.bam

samtools flagstat /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.realign.recal.srt.bam \
  > /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1_final_bam.flagstat

### VARIANT CALLING ############################################################
## Variant calling
echo "Calling variants with GATK HaplotypeCaller..."
/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Xmx12g -jar /usr/local/GenomeAnalysisTK.jar \
  -T HaplotypeCaller \
  -R /PATH/genome.fa \
  -I /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_process_bam/Test1.srt.fixmate.srt.markdup.realign.recal.srt.bam \
  -D /PATH/dbsnp154_hg38.vcf.gz \
  -L /PATH/Exome_RefSeq_targets_hg38.bed
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.snps.indels.raw.without.interval.vcf

### VARIANT FILTERING ##########################################################

echo "=== VARIANT FILTERING PHASE ==="

## Separate SNPs
echo "Separating SNPs..."
java -jar /usr/local/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R /PATH/genome.fa \
  -V /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.snps.indels.raw.without.interval.vcf \
  -selectType SNP \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.raw.snps.vcf

## Separate INDELs
echo "Separating INDELs..."
java -jar /usr/local/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R /PATH/genome.fa \
  -V /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.snps.indels.raw.without.interval.vcf \
  -selectType INDEL \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.raw.indels.vcf

## SNPs filtering
echo "Filtering SNPs..."
java -Xmx12g -jar /usr/local/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R /PATH/genome.fa \
  -V /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.raw.snps.vcf \
  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filterName "my_snp_filter" \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.filtered_snps.vcf

## INDELs filtering
echo "Filtering INDELs..."
java -Xmx12g -jar /usr/local/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R /PATH/genome.fa \
  -V /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.raw.indels.vcf \
  --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
  --filterName "my_indel_filter" \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.filtered_indels.vcf

## Merge filtered VCFs
echo "Merging filtered VCFs..."
java -Xmx12g -jar /usr/local/GenomeAnalysisTK.jar \
  -T CombineVariants \
  -R /PATH/genome.fa \
  --variant:snps /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.filtered_snps.vcf \
  --variant:indels /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.filtered_indels.vcf \
  -o /PATH/exome_analysis_results_hg38/Test1_wes_analysis/Test1_variants_calls/Test1.BOTH.filtered.vcf \
  -genotypeMergeOptions PRIORITIZE \
  -priority foo,bar


### COMPLETION #################################################################

echo "Pipeline completed successfully - $(date)"
