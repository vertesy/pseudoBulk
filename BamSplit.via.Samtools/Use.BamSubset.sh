#!/usr/bin/bash


#Setup
# cd ~/bin
# subl BamSubset.pl
# chmod +x BamSubset.pl
# export PATH=~/bin/:$PATH

#Usage
BAM="/scratch/bioinfo/thomas/Knoblich.GRP/Abel.Vertesy/114593/114593_premRNA/outs/possorted_genome_bam.bam"
CLUSTER="/scratch/bioinfo/thomas/Knoblich.GRP/Abel.Vertesy/114593/114593_premRNA/outs/analysis/clustering/SampleOrigAbel/SampleOrigin.csv"

mkdir BamSplitBySampleOrigin
cd BamSplitBySampleOrigin
date
ml samtools/1.9-foss-2018b
BamSubset.pl -b $BAM -c $CLUSTER
rm Cluster.bam
date
