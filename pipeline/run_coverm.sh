#!/usr/bin/bash
# author: Jinxin Meng
# modified date: 2025-03-15, 22:24:47

( [ $# -lt 3 ] ) && echo -e "Usage: bash $0 [*sort.bam] [genome_filepath] [out_prefix] [covered_fraction: 0 <0-100>]" && exit 127

( [ -f $3.cvm ] ) && echo "Skip sample: $3." && exit 0

cvg=${4:-0}

if [ $cvg -eq 0 ]; then
    coverm genome --bam-files $1 --genome-fasta-list $2 --min-covered-fraction 0 --output-file $3.cvm \
        --methods length count covered_fraction covered_bases tpm rpkm relative_abundance || \
        { rm $3.cvm && exit -1; }
else
    coverm genome --bam-files $1 --genome-fasta-list $2 --min-covered-fraction $4 --output-file $3.cvm \
        --methods covered_fraction covered_bases tpm rpkm relative_abundance || \
        { rm $3.cvm && exit -1; }
fi
