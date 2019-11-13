#!/bin/bash

# usage meta.sh -r <raw> -t <trimmed> -b <bam> -o <out>

echo "Running from $HOSTNAME";

help_msg="Extracts the genes from the specified region and scaffolds,\n
outputs fasta file of the specified gene.\n
Only works for single .fasta files. Output of ~brian/bin/split\n
\n
Usage:\n
    $(basename "$0") [-h] [-r RAW_READS] [-t TRIMMED_READ] [-b BAM_DIR] \n
                     [-o OUTDIR] scaffold_consensus.fasta\n
\n
where:\n
    -h  shows help message\n
    -r  directory containing raw reads\n
    -t  directory containing trimmed reads\n
    -b  directory containing mapped bam files (aligned and coordinate-sorted)\n
    -o  directory to write output files (default: current directory)"

# default output directory lol
OUT=$(pwd)

# getopts does not support optional arguments
# BIG SAD
while getopts ":h:r:b:t:o:" option; do
    case "$option" in
        h) 
            echo -e $help_msg >&2
            exit
            ;;
        r) 
            RAW=$OPTARG >&2
            ;;
        t) 
            TRIM=$OPTARG >&2
            ;;
        b) 
            BAM=$OPTARG >&2
            ;;
        o)
            OUT=$OPTARG >&2
            ;;
        \?)
            echo "Invalid option. -$OPTARG" >&2
            echo -e $help_msg
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            echo -e $help_msg >&2 
            exit 1
            ;;
    esac
done

if [[ ! -d $OUT ]]
then
    mkdir $OUT
fi

# raw read counts
#for sample in `ls -1 $RAW`
#do
#    ID=$(basename -- "$sample")
#    for file in `ls -1 $RAW/$sample`
#    do
#        COUNT=$(zcat $RAW/$sample/$file | wc -l)
#        COUNT=$(echo "scale=0; $COUNT/4" | bc)
#        echo -e "$file\t$COUNT" >> $OUT/${ID}.raw
#    done
#done

# trimmed read counts
for sample in `ls -1 $TRIM`
do
    ID=$(basename -- "$sample")
    for file in `ls -1 $TRIM/$sample/*1P*`
    do
       COUNT=$(zcat $file | wc -l)
       COUNT=$(echo "scale=0; $COUNT/4" | bc)
       fName=$(basename -- "$file")
       echo -e "$fName\t$COUNT" >> $OUT/${ID}.trim
    done
done

# uniquely mapped reads
# lmao I already remove all the unmapped reads
#samtools view accepted_hits.bam | awk '$5==255{print $0}' | wc -l

# what I should do instead, is take all the uniquely mapped reads
# uniqueness of mapping is a probability so I should take reads with
# a mapping quality of above a certain threshold
# 255 is absolutely sure, but the number ranges from 0 to 255
# 10 to 30 is a good range of qualities to threshold from
#for sample in `ls -1 $BAM`
#do
#    ID=$(basename -- "$sample")
#    for bam in `ls $BAM/$ID/*_Aligned*`
#    do
# #       echo $bam
#        samtools view -q 10 -b $bam > $BAM/$ID/${ID}_unique.bam &
#    done
#done
