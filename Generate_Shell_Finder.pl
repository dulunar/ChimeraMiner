#!/usr/bin/perl -w
#########################################################################
# File Name: Run_Finder.pl
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

my $idir = dirname(File::Spec->rel2abs( $0 )); chomp $idir;

my ($list,$out,$ref,$step1,$step2,$len,$help);

GetOptions
(	
	"o=s"=>\$out,
	"i=s"=>\$list,
	"S1:s"=>\$step1,
	"S2:s"=>\$step2,
	"L=i"=>\$len,
	"r=s"=>\$ref,
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
	-r <complete path of genome reference> <default: /home/luna/Desktop/database/homo_bwa/hsa.fa>

INFO

die $usage if (!$list || $help);

$out ||= "runFinder.sh";
$step1 ||= "$idir/Insertion.SRExtract.ReConFastq.pl";
$step2 ||= "$idir/SearchOverlapSEchimera.pl";
$len ||= 30;
$ref ||= "/home/luna/Desktop/database/homo_bwa/hsa.fa";

`samtools faidx $ref` if(! -e "${ref}.fai"); # check if there indexed for reference.

my $refdir = dirname(File::Spec->rel2abs( $ref ));      chomp $refdir;

my $sdir;

open IN,"< $list" || die $!;
while(<IN>){
	chomp;
	my ($samp,$bam) = split /\s+/,$_;
	my $dir = dirname($bam);
	$sdir = $dir;
	
	$out = "runFinder.$samp.sh";

	open OH,"> $sdir/$out" || die $!;
	
	print OH "echo starts\ndate\n";
	print OH "perl $step1 -i $bam -m $samp -d $dir -r $ref && echo first split and find INSERT chimerasjob done\ndate\n";
	print OH "sh $dir/Chr_split/run.aln.sh &> $dir/Chr_split/run.aln.sh.log && echo realign soft-alignment reads to reference\n";
	print OH "split -l 6 $dir/run.$samp.4Search.sh $dir/run.4Search\n";
	print OH "chmod +x $dir/run.4Search*\n";
	print OH "ls $dir/run.4Search* | perl -ne \'chomp;\`nohup sh \$_ &> \$_.log &\`;\'\n";
	
	close OH;

	open OUT,"> $sdir/run.$samp.4Search.sh" || die $!;
	for my $chr(1..22,"X","Y","MT"){
		my $realigngz = "$dir/Chr_split/$samp.chr$chr.sam.gz";
		print OUT "date\nperl $step2 -i $realigngz -o $samp.chr$chr -L $len -d $dir -r $refdir > log.step4.chr$chr\ngzip $dir/chimeras/$samp.chr$chr\ngzip $dir/chimeras/$samp.chr$chr.err\necho chr$chr is done\ndate\n";
	}
	close OUT;
}
close IN;
#close OH;

`chmod +x $sdir/*.sh`;
