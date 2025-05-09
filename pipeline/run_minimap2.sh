#!/usr/bin/bash
# author: Jinxin Meng
# modified date: 2025-03-15, 22:24:47

( [ $# -ne 3 ] ) &&  { echo -e "Usage: $0 [fq|fq1,fq2] [*.mmi|*.fa] [out_prefix]" && exit 2; }
( [ -e $3.sort.bam ] ) && { echo -e "Skip sample: ${3##*/} .." && exit 0; }

if [[ $1 =~ "," ]];then
    fq1=$(echo $1 | cut -d "," -f1)
    fq2=$(echo $1 | cut -d "," -f2)
    minimap2 --MD -t 24 -ax sr $2 $fq1 $fq2 2>>$3.log |\
        samtools view -F 4 -@ 24 -bS | samtools sort -@ 24 -o $3.sort.bam -
else
    minimap2 --MD -t 24 -ax sr $2 $1 2>>$3.log |\
        samtools view -F 4 -@ 24 -bS | samtools sort -@ 24 -o $3.sort.bam -
fi
