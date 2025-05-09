#!/usr/bin/bash
# author: Jinxin Meng
# modified date: 2025-03-15, 22:24:47

if [ $# -ne 3 ];then
    echo -e "$0 [fq|fq1,fq2] [host] [out_prefix]\e[0m"
    echo -e " \033[31mOptional host:\033[0m:"
    printf "%10s:    %-20s %-20s\n" human /data/database/host_genomes/Homo_GRCh38.p14/GCF_000001405.40_GRCh38.p14
    printf "%10s:    %-20s %-20s\n" pig /data/database/host_genomes/Sus_scrofa11.1/GCF_000003025.6_Sscrofa11.1
    printf "%10s:    %-20s %-20s\n" mouse /data/database/host_genomes/Mus_GRCm39/GCF_000001635.27_GRCm39
    exit 2
fi

fq=$1; host=$2; out=$3; trds=16
declare -A host_index
host_index["human"]="/data/database/host_genomes/Homo_GRCh38.p14/GCF_000001405.40_GRCh38.p14"
host_index["mouse"]="/data/database/host_genomes/Mus_GRCm39/GCF_000001635.27_GRCm39"
host_index["pig"]="/data/database/host_genomes/Sus_scrofa11.1/GCF_000003025.6_Sscrofa11.1"
index=${host_index["${host}"]}

( [ -f ${out}_map2host.log ] && grep -q 'overall' ${out}_map2host.log ) && echo -e "Skip sample: $out." && exit 0

if [[ ${fq} =~ "," ]];then
    fq1=$(echo $fq | cut -d "," -f1)
    fq2=$(echo $fq | cut -d "," -f2)
    fastp -w $trds -q 20 -u 30 -n 5 -y -Y 30 -l 80 --trim_poly_g \
        -i $fq1 -I $fq2 -o ${out}_clean_1.fq.gz -O ${out}_clean_2.fq.gz -h /dev/null -j /dev/null 1>${out}_fastp.log 2>&1
    bowtie2 --end-to-end --mm --fast -p $trds -x $index --no-head -1 ${out}_clean_1.fq.gz -2 ${out}_clean_2.fq.gz 2> ${out}_map2host.log |\
        perl -ne 'chomp;@s=split /\s+/;if($s[1]==77){print "\@$s[0]/1\n$s[9]\n+\n$s[10]\n";}elsif($s[1]==141){print STDERR "\@$s[0]/2\n$s[9]\n+\n$s[10]\n";}' \
        > >(pigz> ${out}_rmhost_1.fq.gz) 2> >(pigz > ${out}_rmhost_2.fq.gz)
    chmod 444 ${out}_fastp.log ${out}_map2host.log ${out}_rmhost_1.fq.gz ${out}_rmhost_2.fq.gz
    rm ${out}_clean_1.fq.gz ${out}_clean_2.fq.gz
else
    fastp -w $trds -q 20 -u 30 -n 5 -y -Y 30 -l 80 --trim_poly_g \
        -i $fq -o ${out}_clean.fq.gz -h /dev/null -j /dev/null 1>${out}_fastp.log 2>&1
    bowtie2 --end-to-end --mm --fast -p $trds -x $index --no-head -U ${out}_clean.fq.gz 2> ${out}_map2host.log |\
        perl -ne 'chomp;@s=split /\s+/;if($s[1]==4){print "\@$s[0]\n$s[9]\n+\n$s[10]\n";}' \
        > >(pigz > ${out}_rmhost.fq.gz)
    chmod 444 ${out}_fastp.log ${out}_map2host.log ${out}_rmhost.fq.gz
    rm ${out}_clean.fq.gz
fi
# minimap2 比对的序列更多，保留的序列更少
