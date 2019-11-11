#!/bin/bash

#Usage: mapping.sh sampleDir resultsDir genomeDir
# input arguments: sample directory, name of sample, and directory in which the reads are stored
sampleDir=$1

resultsDir=$2

genomeDir=$3

if [ ! -d "$sampleDir" ]; then
    echo "Samples directory does not exist";
    exit
fi

# 1. 1pass
# generate reference genome
#STAR --runMode genomeGenerate --genomeDir $genomeDir \
# --sjdbGTFfile $genomeDir/*.gtf \
# --genomeFastaFiles $genomeDir/*.fa --runThreadN 5
    
for sample in `ls -1 $sampleDir`; do
    #get names of files for naming!
    ID=$(basename -- "$sample")
    echo $ID
    # make directories for analysis
    mkdir -p $resultsDir/1pass/$ID/ \
             $resultsDir/2pass/$ID/ \
             $resultsDir/2pass/$ID/genome/ \
             $resultsDir/2pass/$ID/outfiles/
    echo "1pass"
    R1=$( ls $sampleDir/$sample/*1P* | tr '\n' ',' | sed 's/.$//' )
    R2=$( ls $sampleDir/$sample/*2P* | tr '\n' ',' | sed 's/.$//' )
#    echo "$R1 spaaaace  $R2"
    # separate paired-end mates by a space
    # separate mutilple read files by a comma 
    # (with no spaces in between)
    STARlong --genomeDir $genomeDir \
             --readFilesIn $R1 $R2 \
             --readFilesCommand zcat \
             --runThreadN 15 \
             --limitBAMsortRAM 20000000000 \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix $resultsDir/1pass/$ID/${ID}_ 2>&1 $resultsDir/2pass/$ID/outfiles/0-1pass.out
    echo "done"
    echo "2pass"
    # 2. 2pass
    echo "Genome gen"
    STAR --runMode genomeGenerate \
         --genomeDir $resultsDir/2pass/$ID/genome/\
         --sjdbGTFfile $genomeDir/*.gtf \
         --genomeFastaFiles $genomeDir/*.fa \
         --sjdbFileChrStartEnd $resultsDir/1pass/$ID/${ID}_SJ.out.tab \
         --sjdbOverhang 75 \
         --runThreadN 5 2>&1 $resultsDir/2pass/$ID/outfiles/1-2pass.out

    echo "2pass map"
    STARlong --genomeDir $resultsDir/2pass/$ID/genome/ \
             --readFilesIn $R1 $R2 \
             --readFilesCommand zcat \
             --runThreadN 15 \
             --limitBAMsortRAM 20000000000 \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix $resultsDir/2pass/$ID/${ID}_ 2>&1 $resultsDir/2pass/$ID/outfiles/2-2passalign.out
    echo "done"
done 

