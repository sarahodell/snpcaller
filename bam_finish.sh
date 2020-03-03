#!/bin/bash
### Pipeline for bam file processing post-bwa mem in order to call SNPs

module load jdk/1.8
module load R
module load java/1.8
module load maven
module load bwa/0.7.17.r1188
module load picardtools/2.7.1
module load GATK/3.6

picard=/share/apps/picard-tools-2.7.1/picard.jar
GATK=/share/apps/GATK-3.6/GenomeAnalysisTK.jar
genome=/home/sodell/projects/resources/Zea_maysAGPv3/Zea_mays.AGPv3.dna.toplevel.fa
sample="$1"

java -Xmx128g -jar $GATK \
     -T IndelRealigner \
     -R $genome \
     -I "$sample"_sorted_markdup.bam \
     -targetIntervals "$sample"_realignment_targets.list \
     -o "$sample"_removedup_realigned.bam

#java -Xmx128g -jar $GATK \
#     -T BaseRecalibrator \
#     -R $genome \
#     -I "$sample"_removedup_realigned.bam \
#     -o "$sample"_removedup_realigned_recal.table

#java -Xmx128g -jar $GATK \
#     -T PrintReads \
#     -R $genome \
#     -I "$sample"_removedup_realigned.bam \
#     -BQSR "$sample"_removedup_realigned_recal.table \
#     -o "$sample"_removedup_religned_recal.bam

#java -Xmx128g -jar $GATK \
#     -T BaseRecalibrator \
#     -R $genome \
#     -I "$sample"_removedup_realigned.bam \
#     -BQSR "$sample"_removedup_realigned_recal.table \
#     -o "$sample"_after_recal.table

#java -Xmx128g -jar $GATK \
#     -T AnalyzeCovariates \
#     -R $genome \
#     -after "$sample"_after_recal.table \
#     -plots "$sample"_recal_plots.pdf

#java -Xmx128g -jar $picard SortSam \
#     INPUT="$sample"_bwa_mem.sam \
#     OUTPUT="$sample".bam \
#     SORT_ORDER=coordinate

#java -Xmx40g -jar $picard CollectAlignmentSummaryMetrics \
#     R=$genome \
#     I="$sample".bam \
#     O="$sample"_alignment_metrics.txt

#java -Xmx40g -jar $picard CollectInsertSizeMetrics \
#     INPUT="$sample".bam \
#     OUTPUT="$sample"_insert_metrics.txt \
#     HISTOGRAM_FILE="$sample"_insert_size_histogram.pdf

#samtools depth -a "$sample".bam > "$sample"_depth_out.txt

#java -Xmx40g -jar $picard AddOrReplaceReadGroups \
#     I="$sample".bam \
#     O="$sample"_with_RG.bam \
#     SORT_ORDER=coordinate \
#     RGID="$sample" \
#     RGLB="$sample" \
#     RGPL="illumina" \
#     RGPU="$sample" \
#     RGSM="$sample"

#java -Xmx40g -jar $picard MarkDuplicates \
#     VALIDATION_STRINGENCY=LENIENT \
#     INPUT="$sample"_with_RG.bam \
#     OUTPUT="$sample"_sorted_markdup.bam \
#     METRICS_FILE="$sample"_duplicate_metrics.txt

#java -Xmx40g -jar $picard BuildBamIndex \
