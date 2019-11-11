#!/bin/bash

#usage trim.sh <reads> <out>

TRIM=/usr/local/trimmomatic

READS=$1

OUT=$2

if [[ ! -d $OUT ]]
then
    mkdir $OUT
fi

for i in `ls $READS`
do 
    if [[ ! -d $OUT/$i ]]
    then
        mkdir $OUT/$i
    fi
    R1=$(ls $READS/$i/*R1*)
    R2=$(ls $READS/$i/*R2*)
    paste <(echo -e "$R1") <(echo -e "$R2") > $OUT/tmp

    while read read1 read2
        do
            NAME=$(basename -- "$read1")
            NAME="${NAME%.fastq.gz}"
            #            echo $NAME
            echo -e "Trimming $read1 and $read2"
            java -jar $TRIM/trimmomatic-0.36.jar PE -phred33 -threads 4 \
                -trimlog $OUT/$i/${NAME}_trim.log $read1 $read2 -baseout $OUT/$i/$NAME.fq.gz\
            ILLUMINACLIP:$TRIM/adapters/TruSeq2-PE.fa:2:30:10:5:TRUE\
                LEADING:3 TRAILING:3 SLIDINGWINDOW:3:30 MINLEN:36
            echo -e "Done"
        done < $OUT/tmp
done
