#!/bin/bash

# This script for snp calling

set -eu
if [ "$#" -ne 6 ]; then
	echo "Usage: SNP_calling.sh <software pathway> <sample name> <reference fasta> <fastq1> <fastq2> <output directory>"
	exit 1
fi

SOFTWARE=$1
SAMPLE=$2
REFERENCE=$3
FQ1=$4
FQ2=$5
OUTDIR=$6

# Pathway of software, including gatk, samtools, bwa, BamDeal
#SOFTWARE="/path/to/software/" 

echo "mkdir -p -m 755 $OUTDIR/$SAMPLE"
echo "mkdir -p -m 755 $OUTDIR/$SAMPLE/java_tmp"
echo "$SOFTWARE/bwa mem -t 5 -M -R \"@RG\tID:$SAMPLE\tPL:illumina\tPU:$SAMPLE\t\tLB:$SAMPLE\tSM:$SAMPLE\tCN:BGI\" $REFERENCE $FQ1 $FQ2 | $SOFTWARE/samtools view -Sb -o $OUTDIR/$SAMPLE/${SAMPLE}.bam -"
echo "$SOFTWARE/java -Xmx2G -jar $SOFTWARE/SortSam.jar INPUT=$OUTDIR/$SAMPLE/${SAMPLE}.bam TMP_DIR=$OUTDIR/$SAMPLE/java_tmp OUTPUT=$OUTDIR/$SAMPLE/${SAMPLE}.sort.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
echo "ln -sf $REFERENCE $OUTDIR/$SAMPLE/${SAMPLE}.sort.fa"
echo "$SOFTWARE/samtools faidx $OUTDIR/$SAMPLE/${SAMPLE}.sort.fa"
echo "$SOFTWARE/java -jar $SOFTWARE/CreateSequenceDictionary.jar R=$OUTDIR/$SAMPLE/${SAMPLE}.sort.fa O=$OUTDIR/$SAMPLE/${SAMPLE}.sort.dict"
echo "$SOFTWARE/samtools index $OUTDIR/$SAMPLE/${SAMPLE}.sort.bam"
echo "$SOFTWARE/java -Xmx3G -Djava.io.tmpdir=$OUTDIR/$SAMPLE/java_tmp -jar $SOFTWARE/MarkDuplicates.jar REMOVE_DUPLICATES=false I=$OUTDIR/$SAMPLE/${SAMPLE}.sort.bam O=$OUTDIR/$SAMPLE/markDup.bam METRICS_FILE=$OUTDIR/$SAMPLE/${SAMPLE}markDup.metrics VALIDATION_STRINGENCY=SILENT"
echo "$SOFTWARE/samtools index $OUTDIR/$SAMPLE/markDup.bam"
echo "$SOFTWARE/BamDeal statistics Coverage -i $OUTDIR/$SAMPLE/${SAMPLE}.markDup.bam -r $OUTDIR/$SAMPLE/${SAMPLE}.sort.fa -o ${SAMPLE}"
echo "$SOFTWARE/java -Xmx3G -Djava.io.tmpdir=$OUTDIR/$SAMPLE/java_tmp -jar $SOFTWARE/GenomeAnalysisTK.jar -T UnifiedGenotyper -nt 12 -R $OUTDIR/$SAMPLE/${SAMPLE}.sort.fa -I $OUTDIR/$SAMPLE/${SAMPLE}.markDup.bam -allowPotentiallyMisencodedQuals -stand_call_conf 20.0 -o $OUTDIR/$SAMPLE/${SAMPLE}.gatk.vcf"
echo "$SOFTWARE/java -Xmx3G -Djava.io.tmpdir=$OUTDIR/$SAMPLE/java_tmp -jar $SOFTWARE/GenomeAnalysisTK.jar -T SelectVariants -R $OUTDIR/$SAMPLE/${SAMPLE}.sort.fa -selectType SNP -V $OUTDIR/$SAMPLE/${SAMPLE}.gatk.vcf -o $OUTDIR/$SAMPLE/${SAMPLE}.gatk.select_snp.vcf"
echo "$SOFTWARE/java -Xmx3G -Djava.io.tmpdir=$OUTDIR/$SAMPLE/java_tmp -jar $SOFTWARE/GenomeAnalysisTK.jar -T VariantFiltration -R $OUTDIR/$SAMPLE/${SAMPLE}.sort.fa -V $OUTDIR/$SAMPLE/${SAMPLE}.gatk.select_snp.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o $OUTDIR/$SAMPLE/${SAMPLE}.raw.snp.vcf"
