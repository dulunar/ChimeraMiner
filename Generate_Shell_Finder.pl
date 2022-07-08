#!/usr/bin/perl -w
#########################################################################
# File Name: Generate_Shell_Finder.pl
# Author: Luna
# Mail: nlu@seu.edu.cn
# Created Time: Tue 25 Sep 2018 06:44:00 PM CST
#########################################################################
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use File::Spec;

our $AUTHOR = '$Author: Na Lu <nlu@seu.edu.cn> $';

my $red = "\033[0;31m";
my $end = "\033[0m";

my $idir = dirname(File::Spec->rel2abs( $0 )); chomp $idir;

my ($list,$out,$ref,$step1,$step2,$len,$samtools, $bwa, $help);

GetOptions
(	
	"o=s"=>\$out,
	"i=s"=>\$list,
	"S1:s"=>\$step1,
	"S2:s"=>\$step2,
	"L=i"=>\$len,
	"r=s"=>\$ref,
	"b=s"=>\$bwa,
	"t=s"=>\$samtools,
	"help|?"=>\$help,
);

my $usage=<<INFO;

Usage:
	perl $0 [options]
Options:

	-i <file> <in.input file, two columns in this file, SampleName	DirectoryBam>
	-o <file> <the work shell file for running finder> <default:runFinder.sh>
	-S1 <file> <the first scripts for Finder> <default:$idir/Insertion.SRExtract.ReConFastq.pl>
	-S2 <file> <the second scripts for Finder> <default:$idir/SearchOverlapSEchimera.pl>
	-L <INT> <min Length of each segment> <default:30>
	-r <${red}complete path of genome reference$end> <default: /home/luna/Desktop/database/homo_bwa/hsa.fa>
	-b <${red}the complete path of bwa> (dflt: `which bwa`)$end
	-t <${red}the complete path of samtools> (dflt: `which samtools`)$end

INFO

die $usage if (!$list || $help);

$out ||= "runFinder.sh";
$step1 ||= "$idir/Insertion.SRExtract.ReConFastq.pl";
$step2 ||= "$idir/SearchOverlapSEchimera.pl";
$len ||= 30;
$ref ||= "/home/luna/Desktop/database/homo_bwa/hsa.fa";

die "$red$ref did not exists$end\n" if(! -e "$ref");

if(! defined $samtools){
	$samtools = `which samtools`;
	$samtools=~s/^s+|\s+$//g;
}
die "$red samtools not exists, please write your samtools path to your environment configuration, maybe ~/.bashrc$end\n" if(! $samtools);
die $red."$samtools not exists, Please Check it$end\n" if(! -e $samtools);

if(! defined $bwa){
	$bwa = `which bwa`;
	$bwa=~s/^\s+|\s+$//g;
}
die "$red bwa did not exists, please write your bwa path to your environment configuration, maybe ~/.bashrc$end\n" if(! $bwa);
die $red."$bwa not exists, Please Check it$end\n" if(! -e $bwa);

`$samtools faidx $ref` if(! -e "${ref}.fai"); # check if there indexed for reference.

my $refdir = dirname(File::Spec->rel2abs( $ref ));
chomp $refdir;

my $sdir;

open IN,"< $list" || die $!;
while(<IN>){
	chomp;
	my ($samp,$bam) = split /\s+/,$_;
	my $dir = dirname($bam);
	$sdir = $dir;

	my $sort = `$samtools view -H $bam | grep "SO:coordinate"`;
	$sort=~s/^\s+|\s+$//g;
	die $red."$bam is a sorted bam, can not used in this pipeline$end" if($sort);

	$out = "runFinder.$samp.sh";

	open OH,"> $sdir/$out" || die $!;
	
	print OH "echo starts\ndate\n";
	print OH "perl $step1 -i $bam -m $samp -d $dir -r $ref -t $samtools -b $bwa && echo first split and find INSERT chimerasjob done\ndate\n";
	print OH "sh $dir/Chr_split/run.aln.sh &> $dir/Chr_split/run.aln.sh.log && echo realign soft-alignment reads to reference\n";
	print OH "split -l 4 $dir/run.$samp.4Search.sh $dir/run.4Search\n";
	print OH "chmod +x $dir/run.4Search*\n";
	print OH "for i in `ls $dir/run.4Search* | grep -v log`\ndo\n\tnohup sh \$i &> \$i.log &\ndone\n";
	
	close OH;

	open OUT,"> $sdir/run.$samp.4Search.sh" || die $!;
	for my $chr(1..22,"X","Y","MT"){
		my $realigngz = "$dir/Chr_split/$samp.chr$chr.sam.gz";
		print OUT "date\nperl $step2 -i $realigngz -o $samp.chr$chr -L $len -d $dir -r $refdir -t $samtools\n";
		print OUT "echo chr$chr is done\ndate\n";
	}
	close OUT;
}
close IN;

`chmod +x $sdir/*.sh`;
