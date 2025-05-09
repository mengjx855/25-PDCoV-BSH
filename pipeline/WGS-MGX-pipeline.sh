# qc
cat raw_fq.filelist | parallel -j 6 --colsep="\t" run_fastp_rmhost.sh {2} pig clean_fq/{1}


# assemble
cat clean_fq.filelist | parallel -j 5 --colsep="\t" run_megahit.sh {2} assembly/{1}
cat ../sample_name | parallel -j 2 seqkit seq -g -m 1500 {}.fa \| seqkit replace -p \" .\*\" -r \"\" \| seqkit replace -p \"\^\" -r \"{}\|\" -o {}.m1500.fa


# binning
cat sample_name | parallel -j 6 run_minimap2.sh clean_fq/{}_rmhost_1.fq.gz,clean_fq/{}_rmhost_2.fq.gz assembly/{}.m1500.fa semibin2/{}
cat sample_name | parallel -j 6 SemiBin2 single_easy_bin -t 32 --random-seed 2025 -i assembly/{}.m1500.fa --input-bam semibin2/{}.sort.bam --environment pig_gut -o semibin2/{}
ls */output_bins/* | perl -ne 'chomp;$x=$_;@s=split/\//;$s[2]=~s/SemiBin_//;print "cp $_ bins/$s[0].$s[2]\n"' | sh

# checkm2
checkm2 predict --input semibin2/bins/ --output-directory checkm2 --tmpdir /tmp -x .fa --threads 56 --force
awk -F '\t' '$2>50 && $3<=10' quality_report.tsv > bins.after_qc.MQ.tsv
awk -F '\t' '$2>90 && $3<=5' quality_report.tsv > bins.after_qc.HQ.tsv
awk -F '\t' '$2>=70 && $3<=5 && ($2-$3*5)>55' quality_report.tsv > bins.after_qc.QS.tsv


# dRep
dRep compare dRep/ -g checkm2/bins.after_qc.MQ.filepath -d -pa 0.9 -sa 0.95 -nc 0.30 -cm larger -p 50 --S_algorithm fastANI
parse_dRep.pl dRep/data_tables/Cdb.csv dRep.clu.info -ckm2 checkm2/bins.after_qc.MQ.tsv


# taxa
gtdbtk classify_wf --cpus 80 --pplacer_cpus 12 --skip_ani_screen --genome_dir species_fa/ --out_dir gtdbtk_out -x fa


# abundance
cat ../gtdbtk/species.filelist | parallel -j 5 --colsep="\t" seqkit replace -p \"\.\*\" -r \"{1}\|{nr}\" {2} -o species_fa/{1}.fa
minimap2 -d species.mmi species.fa
cat ../clean_fq.filelist | parallel -j 6 --colsep="\t" ../run_minimap2.sh {2} species.mmi sort_bam/{1}
cat ../sample_name | parallel -j 6 ../run_coverm.sh sort_bam/{}.sort.bam species_fa.filepath cvm_out_c50/{}.c50 50
cat ../sample_name | parallel -j 3 ../run_coverm.sh sort_bam/{}.sort.bam species_fa.filepath cvm_out/{}
combine_file_zy_folder_allsample.py -D cvm_out -suffix .cvm --skip 2 -n 1 -v 3 -o species.rc
combine_file_zy_folder_allsample.py -D cvm_out -suffix .cvm --skip 2 -n 1 -v 4 -o species.cvg
combine_file_zy_folder_allsample.py -D cvm_out -suffix .cvm --skip 2 -n 1 -v 6 -o species.tpm
combine_file_zy_folder_allsample.py -D cvm_out -suffix .cvm --skip 2 -n 1 -v 8 -o species.rela_ab


# mpa4
cat clean_fq.filelist | parallel -j 5 --colsep="\t" run_metaphlan.sh {2} mpa4/{1}
combine_file_zy_folder_allsample.py -D mpa4/ -suffix .tsv --skip 5 -n 1 -v 3 -o mpa4.tsv
parse_mpa.py -i mpa4.tsv -o mpa4 -l p,g,s --short


# prokka
# designate PDCoV_D1_5751.65,PDCoV_D3_5751.132,PDCoV_D3_5753.72,PDCoV_D3_5756.119 as archaea
cat checkm2/bins.after_qc.MQ.filepath | parallel -j 10 echo prokka {} --outdir prokka/{/.} --prefix {/.} --kingdom Bacteria --metagenome --addgenes --cpus 8 --quiet > run_prokka.sh


# prodigal
cat ../checkm2/bins.after_qc.MQ.filelist | parallel -j 5 --colsep="\t" seqkit replace -p \"\.\*\" -r \"{1}\|{nr}\" {2} -o fna/{1}.fa
realpath fna/* | parallel -j 30 prodigal -q -p single -f gff -i {} -o gff/{/.}.gff -a faa/{/.}.faa -d ffn/{/.}.ffn


# geneset
cat ../prodigal/ffn/*ffn | seqkit seq -g -m 100 | seqkit replace -p "\s.*" -r "" > genes.ffn
mmseqs easy-cluster genes_total.ffn genes . --cluster-mode 2 --cov-mode 1 --min-seq-id 0.95 -c 0.9 --kmer-per-seq-scale 0.8 --threads 112


# KEGG diamond
diamond blastp -d /data/database/KEGG/KEGG20230401.dmnd --outfmt 6 --min-score 60 --query-cover 50 --max-target-seqs 10 -p 112 -q mmseqs/geneset.faa -o kofam2/kegg.btp >/dev/null 2>/dev/null
perl -ne 'chomp;@s=split/\t/;if($s[0] ne $a){$s[1]=~/\|(.*)/;print "$s[0]\t$1\n";$a=$s[0]}' ko.btp > ko.tsv
profile_aggregate_kegg.pl kofam2/ko.tsv geneset/geneset.tpm kofam2/diamond.ko.tpm

cat ../prodigal/faa.filelist | parallel -j 20 --colsep="\t" -k diamond blastp -d /data/database/KEGG/KEGG20230401.dmnd --outfmt 6 --min-score 60 --query-cover 50 --max-target-seqs 10 -p 6 -q {2} -o diamond_MAGs/{1}.btp \>/dev/null 2\>/dev/null &
seqkit replace -X ../prodigal/faa.filepath -p " .*" -r "" > MAGs.faa
diamond blastp -d /data/database/KEGG/KEGG20230401.dmnd --outfmt 6 --min-score 60 --query-cover 50 --max-target-seqs 10 -p 112 -q MAGs.faa -o MAGs.btp >/dev/null 2>/dev/null &
grep 'K01442\|K15870' MAGs.tsv | perl -ne 'chomp;$_=~/(\S+)\|/;print "$_\t$1\n"' | csvtk join --left-join -t -H -f "3;2" - ../dRep.clu.spread | csvtk join --left-join -t -H -f "4;1" - ../gtdbtk.tsv > bsh.tsv


# geneset abundance
minimap2 -d geneset.mmi ../mmseqs/geneset.ffn
cat ../clean_fq.filelist | parallel -j 6 --colsep="\t" ../run_minimap2.sh {2} geneset.mmi sort_bam/{1}
cat ../sample_name | parallel -j 5 samtools coverage sort_bam/{}.sort.bam -o cvg/{}.cvg
combine_file_zy_folder_allsample.py -D cvg/ -suffix .cvg -n 1 -v 4 -t 1 -o geneset.rc
combine_file_zy_folder_allsample.py -D cvg/ -suffix .cvg -n 1 -v 6 -t 1 -o geneset.cvg
# profile_filter_by_cvg.py geneset.rc geneset.cvg 10 geneset.rc.f
profile_rc2tpm.pl geneset.rc ../mmseqs/geneset.len geneset.tpm


# phylophlan
phylophlan_write_config_file -d a -o phylophlan/phylophlan.cfg --db_aa diamond --map_aa diamond --msa mafft --trim trimal --tree1 fasttree --overwrite
phylophlan -i prodigal/faa_species -t a -o phylophlan/ -d phylophlan --databases_folder /data/database/phylophlan -f phylophlan/phylophlan.cfg --diversity high --nproc 100

