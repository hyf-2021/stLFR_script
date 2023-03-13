#!/bin/bash

# This script for 

set -eu
if [ "$#" -ne 4 ]; then
	echo "Usage: Bins_Dereplication_annotation_gtdbtk.sh <software pathway> <Bins diectory> <output directory> <Bins information file>"
	exit 1
fi

SOFTWARE=$1
BINS=$2
OUTDIR=$3
GENOMEINFO=$4

#source /zfssz5/BC_PS/huangyf/Bin/metagenomics/modules/Env.source.sh
echo "mkdir -m 755 -p $OUTDIR"
echo "mkdir -m 755 -p $OUTDIR/Dereplication_bins"

echo "$SOFTWARE/python $SOFTWARE/dRep dereplicate $OUTDIR -pa 0.9 -sa 0.95 -nc 0.6 -cm larger --genomeInfo $GENOMEINFO"
echo "$SOFTWARE/gtdbtk classify_wf --genome_dir $OUTDIR/Dereplication_bins --out_dir $OUTDIR -x fa --min_perc_aa 10 --cpus 12 --force --pplacer_cpus 12"
