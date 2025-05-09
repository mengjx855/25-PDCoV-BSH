#!/usr/bin/perl
# Jinxin Meng, jinxmeng@zju.edu.cn
# created date: 2023-10-02, 01:49:47
# modified date: 2025-03-06, 08:32:39
use warnings;
use strict;
use Getopt::Long;

my $use = "Usage: perl $0 [OPTIONS] [Cdb.csv] [out_file]
OPTIONS:
  -ckm  : checkm quality socre;
  -ckm2 : checkm2 quality report; 
  -len  : genome length, tab-delimited file with field name and length;
EXAMPLE:
  perl parse_dRep.pl Cdb.csv genome_clu.tsv
  perl parse_dRep.pl -ckm2 quality_report.tsv Cdb.csv genome_clu.tsv\n";

my ($ckm, $ckm2, $len) = (undef, undef, undef);
&GetOptions(
    "ckm:s"  => \$ckm,
    "ckm2:s" => \$ckm2,
    "len:s"  => \$len);

die "$use" unless @ARGV eq 2;
my ($cdb, $out) = @ARGV;

open I, "<$cdb";
my (%clu, @name) = (); 
while (<I>) {
    chomp;
    next if (/primary_cluster/);
    my @s = split /,/;
    $s[0] =~ s/.(fa|fna)//;
    push @{$clu{$s[1]}}, $s[0];
    push @name, $s[0];
}

open O, ">$out";

if ( !defined $ckm && !defined $ckm2 && !defined $len ) {
    for my $i (sort keys %clu) {
        my $count = scalar @{$clu{$i}};
        print O "clu.$i\t$count\t".join(",", @{$clu{$i}})."\n";
    }
    close O;
    exit 0;
}

my %value; 
if (defined $ckm) {
    open I, "$ckm";
    while (<I>) {
        chomp;
        next if /Completeness/;
        my @s = split /\t/;
        $value{$s[0]} = $s[11] - 5 * $s[12] if (grep {$_ eq $s[0]} @name);
    }
} 

if (defined $ckm2) {
    open I,"<$ckm2";
    while (<I>) {
        chomp;
        next if /Completeness/;
        my @s = split /\t/;
        $value{$s[0]} = $s[1] - 5 * $s[2] if (grep {$_ eq $s[0]} @name);
    }
} 

if (defined $len) {
    open I, "<$len";
    while (<I>) {
        chomp;
        next if /len/;
        my @s = split/\t/;
        $value{$s[0]} = $s[1];
    }
} 

select_rep(\%clu, \%value); # 调用子程序并传递哈希的引用

close O; 

sub select_rep {
    my ($clu, $value) = @_;  # 获取哈希的引用
    my %clu = %{$clu};   # 将哈希引用解引用为哈希
    my %value = %{$value};       # 将哈希引用解引用为哈希
    for my $i (sort keys %clu) {
        my $count = scalar @{$clu{$i}};
        if ($count == 1) {
            print O "clu.$i\t1\t@{$clu{$i}}\t@{$clu{$i}}\n";
        } else {
            my $max = -10000;
            my $rep = "";
            for my $j (@{$clu{$i}}) {
                if ($value{$j} >= $max) {
                    $max = $value{$j};
                    $rep = $j;
                }
            }
            print O "clu.$i\t$count\t$rep\t".join(",", @{$clu{$i}})."\n";
        }   
    }
}
