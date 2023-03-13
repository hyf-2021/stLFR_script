#!/bin/bash

# This script is for assembling stLFR data using MetaTrass

set -eu
if [ "$#" -ne 4 ]; then
	echo "Usage: Assembly_MetaTrass.sh <fastq1> <fastq2> <sample prefix> <output directory>"
	exit 1
fi

FQ1=$1
FQ2=$2
SAMPLE=$3
OUTDIR=$4

PYTHON="/path/to/python3"
TRASS="/path/to/MetaTrass/Trass.py"
REF_KRAKEN_DB="/path/to/uhgg_kraken2-db/"
REF_KRAKEN_FA="/path/to/uhgg_kraken2-fa/"
REF_GENOME_INFO="/path/to//all_single_species_genome_size.uhgg.txt"

echo "mkdir -p -m 755 $OUTDIR/$SAMPLE"
echo $PYTHON $TRASS GC -rawfq1 $FQ1 -rawfq2 $FQ2 -outdir $OUTDIR/$SAMPLE -runnow yes
echo $PYTHON $TRASS TB -cleanfq1 $OUTDIR/$SAMPLE/dir1_cleandata/split_reads.1.fq.gz.clean.gz -cleanfq2 $OUTDIR/$SAMPLE/dir1_cleandata/split_reads.2.fq.gz.clean.gz -thread 10 -sample $SAMPLE -ref_db $REF_KRAKEN_DB -genome_size $REF_GENOME_INFO -outdir $OUTDIR/$SAMPLE -runnow yes
echo $PYTHON $TRASS AP -outdir $OUTDIR/$SAMPLE -ref_fa $REF_KRAKEN_FA -thread 10 -parallel 10 -runnow yes
