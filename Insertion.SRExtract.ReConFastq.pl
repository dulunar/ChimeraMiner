#!/usr/bin/perl -w
#########################################################################
# FileName: Insertion.SRExtract.ReConFastq.pl
# Version: 6cd6a087-188a-4094-ace0-ce63b2af80de
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Tue 04 Aug 2019 05:35:01 PM CST
#########################################################################

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Spec;
use File::Basename;

my $red = "\033[0;31m";
my $rend = "\033[0m";

my $usage = "
for first alignment file:
	Step 1:
		seleted first-type chimeric reads	(++ or --).
	Step 2:
		Extract the sequences which are soft-aligned with xxSxxM , xxMxxS or xxSxxMxxS in alignment format from BWA-aligned SAM file.
		Any Soft-aligned fragment mustn't shorter than log4(Chromosome_Length) nt.
	Step 3:
		Reconstruct FASTQ file from all qualified soft-aligned reads. 
		Qualified reads would be cut into two or three part according to the soft-aligned format.
		All reads would be separated into different chromosomes according to the primary alignment results by BWA.
		Cut reads generated in here would be aligned to Hg19 reference by using BWA-sampe mode.

Usage: perl $0 
		-i <bam file>
		-d <output directory>
		-r <complete path of genome reference> <default: /home/luna/Desktop/database/homo_bwa/hsa.fa>
		-g <genome chromosome(optional)>
		-L <chromosome length(optional)>
			${red}options -g -L have to exist together!$rend
		-m <sample name or ID>
		-f <INT> <the max fragment length of a paired-end read> <default:1000>
";

my ($in,$dir,$ref,$genome,$length,$name,$limitdis,$samtools,$bwa,$help);

GetOptions(
	'i=s'	=>	\$in,
	'd=s'	=>	\$dir,
	'r=s'	=>	\$ref,
	'm=s'	=>	\$name,
	'g=s'   =>  \$genome,
	'L=s'   =>	\$length,
	'f=s'	=>	\$limitdis,
	'b=s'	=>	\$bwa,
	't=s'	=>	\$samtools,
	'h|?'	=>	\$help,
);

die "$usage\n" if ($help || !$in ||  !$dir || !$name);

my $count = 0;

$ref ||= "/home/luna/Desktop/database/homo_bwa/hsa.fa";
die "$red$ref did not exists$rend\n" if(! -e "$ref");

`samtools faidx $ref` if(! -e "${ref}.fai"); # check if there indexed for reference.

my $refdir = dirname(File::Spec->rel2abs( $ref ));
chomp $refdir;
$refdir=~s/^\s+|\s+$//g;

if(! defined $samtools){
    $samtools = `which samtools`;
    $samtools=~s/^\s+|\s+$//g;
}
die "$red samtools not exists, please write your samtools path to your environment$rend\n" if(! $samtools);
die $red."$samtools not exists, Please Check it$rend\n" if(! -e $samtools);

if(! defined $bwa){
    $bwa = `which bwa`;
    $bwa=~s/^\s+|\s+$//g;
}
die "$red bwa did not exists, please write your bwa path to your environment$rend\n" if(! $bwa);
die $red."$bwa not exists, Please Check it$rend\n" if(! -e $bwa);

$limitdis ||= 1000;
my (@chr,@genome);

if(!$length || !$genome){
	@genome = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrMT");
	@chr = (0,249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566,16571);
}
else{
	@chr = split /,/,$length;
	@genome = split /,/,$genome;
}

`mkdir -p $dir` if( ! -d $dir);
if($in =~ /\.bam$/){
	`samtools view -H $in | gzip > $dir/$name.first.sam.gz`;
	`samtools view -H $in > header.sam`;
}
elsif($in =~ /\.sam$/){
	`grep "^@" $in | gzip > $dir/$name.first.sam.gz`;
	`grep "^@" $in > $dir/header.sam`;
}
elsif($in =~ /\.sam\.gz$/){
	`zgrep "^@" $in | gzip > $dir/$name.first.sam.gz`;
	`zgrep "^@" $in > $dir/header.sam`;
}
else{
	die $red."this file format is not suit for this tools$rend\n";
}

if($in =~ /\.bam$/){
	open IN, "$samtools view $in | awk '\$6 !~ /H/' | " or die "$!";
}
elsif($in =~ /\.sam$/){
	open IN, "cat $in | awk '\$6 !~ /H/' | " or die "$!";
}
elsif($in =~ /\.sam\.gz$/){
	open IN, "gzip -dc $in |" or die "$!";	
}
else{
	die $red."this file format is not suit for this tools$rend\n";
}

open OA, "| gzip >  $dir/$name.first.chi.gz" || die $!;
open OC, "| gzip >> $dir/$name.first.sam.gz" || die $!;
open OW, "| gzip > $dir/$name.wasted.gz"  || die $!;
open BAM,"| $samtools view -Sb --reference $ref -o $dir/$name.PE.mappable.bam - " || die $!; # for samtools 1.8
# or open BAM, " | sambamba view -S -f bam -T $ref -o $dir/$name.PE.mappable.bam /dev/stdin" || die $!; # for sambamba 0.7.1
# or open BAM, " | sambamba view -S -f bam -o $dir/$name.PE.mappable.bam /dev/stdin" || die $!;

`rm -rf $dir/Chr_split && mkdir -p $dir/Chr_split` if((-d "$dir/Chr_split"));
`mkdir -p $dir/Chr_split` if(!(-d "$dir/Chr_split"));

my $out1;
my $out2;
my $ot1;
my $ot2;
my $OUT1;
my $OUT2;
my %handle = ();

for my $i(1..22,"X","Y","MT"){
	$out1 = "$dir/Chr_split/chr$i.1.fastq";
	$OUT1 = "Fchr".$i;  ##First chr"$i" file
	$out2 = "$dir/Chr_split/chr$i.2.fastq";
	$OUT2 = "Schr".$i;      ##Second chr"$i" file
	open ( $handle{$OUT1}, "| gzip > $out1.gz" ) or die "$!";
	open ( $handle{$OUT2}, "| gzip > $out2.gz" ) or die "$!";
}

print "###$name start\n";

while (<IN>) {
	chomp;
	next if($_ =~ /^@/);
	$count ++;
	print "$count is done\n" if($count % 10000000 == 0);
	my $line = <IN>;
	chomp $line;
	my ($id_1, $flag_1, $chr_1, $pos_1, $format_1, $seq_1, $qua_1) = (split /\t/, $_)[0, 1, 2, 3, 5, 9, 10];
	my ($id_2, $flag_2, $chr_2, $pos_2, $format_2, $seq_2, $qua_2) = (split /\t/, $line)[0, 1, 2, 3, 5, 9, 10];
	if($format_1 =~ /S/ || $format_2 =~ /S/ || $format_1 eq "*" || $format_2 eq "*"){
		&LengthOfSegment($_) if $format_1 =~ /S/;
		&LengthOfSegment($line) if $format_2 =~ /S/;
		print OW "$_\n" if $format_1 eq "*";
		print OW "$line\n" if $format_2 eq "*";
	}

	else{
		my $str_1 = &strand_code($flag_1);
		my $str_2 = &strand_code($flag_2);

		chomp $seq_1;chomp $seq_2;
		my $map1 = length $seq_1;
		my $map2 = length $seq_2;
		my $len1 = ($format_1 =~ /(\d+)D/)?($map1+$1):($map1); #$len1 = ($format_1 =~ /(\d+)I/)?($len1+$1):($len1);
		my $len2 = ($format_2 =~ /(\d+)D/)?($map2+$1):($map2); #$len2 = ($format_2 =~ /(\d+)I/)?($len2+$1):($len2);
		my $end1 = $pos_1 + $len1;  my $end2 = $pos_2 + $len2;
		if($chr_1 ne $chr_2){
			print OC "$_\n$line\n";
			print OA "$id_1\t$chr_1\t$str_1\t$pos_1\t$end1\t$format_1\t$seq_1\t$chr_2\t$str_2\t$pos_2\t$end2\t$format_2\t$seq_2\tdifferentchromosome\n";
		}
		else{
			if($str_1 eq $str_2){
				print OC "$_\n$line\n";
				my $dis = abs($pos_1 - $end2);			
				print OA "$id_1\t$chr_1\t$str_1\t$pos_1\t$end1\t$format_1\t$seq_1\t$chr_2\t$str_2\t$pos_2\t$end2\t$format_2\t$seq_2\tsamestrand\t$dis\n";
			}
			else{
				my $dis = abs($pos_2 + $len2 - $pos_1);
				if(($str_1 eq "+" && $pos_1 > $pos_2) || ($str_1 eq "-" && $pos_1 < $pos_2)){
					print OC "$_\n$line\n";
					print OA "$id_1\t$chr_1\t$str_1\t$pos_1\t$end1\t$format_1\t$seq_1\t$chr_2\t$str_2\t$pos_2\t$end2\t$format_2\t$seq_2\tdiffstrandRF\t$dis\n";
				}
				else{
					if($dis > $limitdis){
						print OC "$_\n$line\n";
						print OA "$id_1\t$chr_1\t$str_1\t$pos_1\t$end1\t$format_1\t$seq_1\t$chr_2\t$str_2\t$pos_2\t$end2\t$format_2\t$seq_2\tdistanceFault\t$dis\n";
					}
					else{
						print BAM "$_\n$line\n";
					}
				}
			}
		}
	}
}
close IN;
close OA;
close OC;
close $OUT1;
close $OUT2;
close BAM;
close OW;

`samtools reheader -P -i header.sam $dir/$name.PE.mappable.bam`;

print "$count is done\n";
print "###distinguish $name end\n";

open SHELL,"> $dir/Chr_split/run.aln.sh" || die $!;

`rm $dir/Chr_split/$name.aln.sam.gz` if(-e "$dir/Chr_split/$name.aln.sam.gz");
print SHELL "echo Start Realignment\ndate\n";
for my $fq1(`ls $dir/Chr_split/chr*.1.fastq*`){
	chomp $fq1;
	my $fq2=$fq1;
	$fq2 =~s/\.1\.fastq/\.2\.fastq/g;
	my ($chr)=(split /\//,$fq1)[-1]=~/(chr.*)\.1\.fastq/;
	print SHELL "bwa aln -t 4 -l 14 $refdir/$chr.fa $fq1 -f $dir/Chr_split/$chr.1.sai\n";
	print SHELL "bwa aln -t 4 -l 14 $refdir/$chr.fa $fq2 -f $dir/Chr_split/$chr.2.sai\n";
	print SHELL "bwa sampe $refdir/$chr.fa $dir/Chr_split/$chr.1.sai $dir/Chr_split/$chr.2.sai $fq1 $fq2 | gzip > $dir/Chr_split/$name.$chr.sam.gz\n";
	print SHELL "rm $dir/Chr_split/$chr.1.sai $dir/Chr_split/$chr.2.sai\n";
}
print SHELL "echo Finished Realignment\ndate\n";
close SHELL;
`chmod +x $dir/Chr_split/run.aln.sh`;
###`sh $dir/Chr_split/run.aln.sh &> $dir/Chr_split/run.aln.sh.log`;

sub Chr{
	my ($chr) = $_[0];
	my $i = -1;
	my $chr_tag;
	foreach my $chr_tmp(@genome){
		$i++;
		$chr_tag = $i if($chr_tmp eq $chr);
	}
	return $chr_tag;
}

sub LengthOfSegment{
	my $line = $_[0];
	my @temp = split /\t/, $line;
	my $chr = $temp[2];
	my $chr_tag = &Chr($chr) + 1;
	if (($temp[5] =~ /S/) && (!($temp[5] =~ /D/)) && (!($temp[5] =~ /I/))) {
		my $min_soft_len = int ( (log($chr[$chr_tag]))/ (log(4)) ) + 1;
		$min_soft_len = 8 if $min_soft_len < 8;
		if ($temp[5] =~ /^(\d+)S(\d+)M$/) {
			my $soft = $1;
			my $map = $2;
			&ReconstructFastq($line) if $soft >= $min_soft_len;
		}
		elsif ($temp[5] =~ /^(\d+)M(\d+)S$/) {
			my $map = $1;
			my $soft = $2;
			&ReconstructFastq($line) if $soft >= $min_soft_len;
		}
		elsif ($temp[5] =~ /^(\d+)S(\d+)M(\d+)S$/) {
			my $soft_1 = $1;
			my $map = $2;
			my $soft_2 = $3;
			&ReconstructFastq($line) if (($soft_1 >= $min_soft_len) && ($soft_2 >= $min_soft_len));
		}
	}
}

sub strand_code{
	my $flag = $_[0];#数组是@_，但是单个参数是$_[0]开始；
	chomp $flag;
	#print "$flag\t";
	$flag = reverse( sprintf("%b",$flag) );
	my @tmp = split //,$flag;
	my $strand = ($tmp[4] == 0)?("+"):("-");
	return $strand;
}

sub seq_trans {
	my $seq = $_[0];
	chomp $seq;
	my $reverse_seq;
	$seq = uc($seq);
	$seq=~tr/ATCG/TAGC/;
	$reverse_seq = reverse $seq;
	return $reverse_seq;
}

sub ReconstructFastq{
	my $line = $_[0];
	my @temp = split /\t/, $line;
	my $id = $temp[0];
	my $flag = $temp[1];
	my $format = $temp[5];
	my $seq = $temp[9];
	my $qua = $temp[10];

	my $str = &strand_code($flag);
	# open continuous output handles for each chromosome in order for the FASTQ writing
	my $chr_tag = $1 if $temp[2] =~ /(chr.*)/;
	$ot1 = "F".$chr_tag;    ##First chr"$i" file
	$ot2 = "S".$chr_tag;    ##Second chr"$i" file

	# FASTQ output
	if($str eq "+"){
		if($temp[5] =~ /^(\d+)S(\d+)M$/){
			my $soft_len = $1;
			my $map_len = $2;
			my $soft_seq = substr($seq, 0, $soft_len);
			my $soft_qua = substr($qua, 0, $soft_len);
			my $map_seq = substr($seq, $soft_len, );
			my $map_qua = substr($qua, $soft_len, );
			print {$handle{$ot1}} "@"."${id}.S 1\n$soft_seq\n+\n$soft_qua\n";
			print {$handle{$ot2}} "@"."${id}.S 2\n$map_seq\n+\n$map_qua\n";
		}
		elsif($temp[5] =~ /^(\d+)M(\d+)S$/) {
			my $map_len = $1;
			my $soft_len = $2;
			my $map_seq = substr($seq, 0, $map_len);
			my $map_qua = substr($qua, 0, $map_len);
			my $soft_seq = substr($seq, $map_len, );
			my $soft_qua = substr($qua, $map_len, );
			print {$handle{$ot2}} "@"."${id}.M 2\n$soft_seq\n+\n$soft_qua\n";
			print {$handle{$ot1}} "@"."${id}.M 1\n$map_seq\n+\n$map_qua\n";
		}
		elsif($temp[5] =~ /^(\d+)S(\d+)M(\d+)S$/) {
			my $soft_len_1 = $1;
			my $map_len = $2;
			my $soft_len_2 = $3;
			my $soft_seq_1 = substr($seq, 0, $soft_len_1);
			my $soft_qua_1 = substr($qua, 0, $soft_len_1);
			my $map_seq = substr($seq, $soft_len_1, $map_len);
			my $map_qua = substr($qua, $soft_len_1, $map_len);
			my $soft_seq_2 = substr($seq, $soft_len_1 + $map_len);
			my $soft_qua_2 = substr($qua, $soft_len_1 + $map_len);
			print {$handle{$ot1}} "@"."$id.1 1\n$soft_seq_1\n+\n$soft_qua_1\n";
			print {$handle{$ot2}} "@"."$id.1 2\n$map_seq\n+\n$map_qua\n";
			print {$handle{$ot1}} "@"."$id.2 1\n$map_seq\n+\n$map_qua\n";
			print {$handle{$ot2}} "@"."$id.2 2\n$soft_seq_2\n+\n$soft_qua_2\n";
		}
	}

	elsif($str eq "-"){
		$seq = &seq_trans($seq); $qua = reverse $qua;
		if($temp[5] =~ /^(\d+)S(\d+)M$/){
			my $soft_len = $1;
			my $map_len = $2;
			my $soft_seq = substr($seq, $map_len);
			my $soft_qua = substr($qua, $map_len);
			my $map_seq = substr($seq, 0, $map_len);
			my $map_qua = substr($qua, 0, $map_len);
			print {$handle{$ot2}} "@"."$id.M 2\n$soft_seq\n+\n$soft_qua\n";
			print {$handle{$ot1}} "@"."$id.M 1\n$map_seq\n+\n$map_qua\n";
		}
		elsif($temp[5] =~ /^(\d+)M(\d+)S$/) {
			my $map_len = $1;
			my $soft_len = $2;
			my $map_seq = substr($seq, $soft_len);
			my $map_qua = substr($qua, $soft_len);
			my $soft_seq = substr($seq, 0, $soft_len);
			my $soft_qua = substr($qua, 0, $soft_len);
			print {$handle{$ot1}} "@"."$id.S 1\n$soft_seq\n+\n$soft_qua\n";
			print {$handle{$ot2}} "@"."$id.S 2\n$map_seq\n+\n$map_qua\n";
		}
		elsif($temp[5] =~ /^(\d+)S(\d+)M(\d+)S$/) {
			my $soft_len_1 = $1;
			my $map_len = $2;
			my $soft_len_2 = $3;
			my $soft_seq_1 = substr($seq, $soft_len_2 + $map_len);
			my $soft_qua_1 = substr($qua, $soft_len_2 + $map_len);
			my $map_seq = substr($seq, $soft_len_2, $map_len);
			my $map_qua = substr($qua, $soft_len_2, $map_len);
			my $soft_seq_2 = substr($seq, 0, $soft_len_2);
			my $soft_qua_2 = substr($qua, 0, $soft_len_2);
			print {$handle{$ot2}} "@"."$id.2 2\n$soft_seq_1\n+\n$soft_qua_1\n";
			print {$handle{$ot1}} "@"."$id.2 1\n$map_seq\n+\n$map_qua\n";
			print {$handle{$ot2}} "@"."$id.1 2\n$map_seq\n+\n$map_qua\n";
			print {$handle{$ot1}} "@"."$id.1 1\n$soft_seq_2\n+\n$soft_qua_2\n";
		}
	}
}
