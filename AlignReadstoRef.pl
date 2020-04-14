#!/usr/bin/perl -w
#########################################################################
# File Name: AlignReadstoRef.pl
# Author: Luna
# Mail: nlu@seu.edu.cn
# Created Time: Fri 21 Oct 2016 21:18:27 CST
#########################################################################
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
my ($fq,$ref,$hdir,$sh,$algo,$dh,$cpu,$help);

GetOptions
(
	"fq=s"=>\$fq,
	"ref=s"=>\$ref,
	"hdir=s"=>\$hdir,
	"Algo=s"=>\$algo,
	"sh=s"=>\$sh,
	"dh=s"=>\$dh,
	"T|cpu=s"=>\$cpu,
	"help|?"=>\$help,
);

my $usage=<<INFO;

Usage:
		perl $0 [options]
		perl $0 -fq /home/luna/work/Chimeras/underwrite/sc_mda/fastq.scMDA2.lst -hdir /home/luna/work/Chimeras/underwrite/sc_mda -sh a01.runbwa -dh yes|no
Options:
		-fq	<STRING>        The fq|fastq file list (Name	FQ1	FQ2)
		-ref    <STRING>        ref.fasta <Default:/home/luna/Desktop/database/homo_bwa/hsa.fa>
		-hdir    <STRING>        home output directory
		-Algo   <STRING>        three alignment algorithms:<Default:mem>|aln/samse
		-sh     <STRING>        output work.sh (Default: a01.run.bwa)
		-cpu/-T	<INT>	The number threads used for aligning (Default: 8)
		-dh	<STRING>	Whether delete Hard clipped Alignment:<Default:no>(yes|no)

INFO

die $usage if ($help || !$fq || !$hdir);
my $bwa = "/home/luna/Desktop/Software/bwa/bwa";
$algo ||= "mem";
$ref ||= "/home/luna/Desktop/database/homo_bwa/hsa.fa";
$sh ||= "a01.run.bwa";
$cpu ||= 20;
open LIST, " $fq  " || die $!;
while(<LIST>){
	
	chomp;
	my ($name,$fq1,$fq2) = (split /\t/,$_)[0,1,2];
	
	my $outdir = "$hdir/$name";
	$outdir = abs_path($outdir);
	`mkdir -p $outdir && chmod 755 $outdir` unless (-e $outdir);
	open OUT, " > $outdir/$sh.$name.sh " || die $!; 
	open BAMLST, " > $outdir/bam.lst " || die $!;
	
	print OUT "echo starts\ndate\n";
	if($dh !~ /yes/){
		print OUT "$bwa mem -t $cpu -k 30 -R '\@RG\\tID:$name\\tLB:$name\\tSM:$name\\tPL:ILLUMINA' $ref $fq1 $fq2 | samtools view -Sb -t $ref.fai -o $outdir/$name.mem.bam - && echo bwa $name job done\n";
		print OUT "/home/luna/Desktop/Software/bamtools/bin/bamtools filter -in $outdir/$name.mem.bam -mapQuality \">=1\" -out $outdir/$name.filter.bam\n";
		print OUT "java -Xmx24g -jar /home/luna/Desktop/Software/picard-tools-1.119/SortSam.jar TMP_DIR=/home/luna/Desktop/work/temp_dir VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=2000000 I=$outdir/$name.filter.bam O=$outdir/$name.sort.bam\n";
		print OUT "java -Xmx15g -jar /home/luna/Desktop/Software/picard-tools-1.119/MarkDuplicates.jar TMP_DIR=/home/luna/Desktop/work/temp_dir VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1020 I=$outdir/$name.sort.bam O=$outdir/$name.final.bam METRICS_FILE=$outdir/$name.filter.dedup.metrics\n";
		print BAMLST "$name\t$outdir/$name.mem.bam\n";
		print BAMLST "$name.filter\t$outdir/$name.filter.bam\n";
		print BAMLST "$name.sort\t$outdir/$name.sort.bam\n";
		print BAMLST "$name.final\t$outdir/$name.final.bam\n";
	}
	else{
		print OUT "$bwa mem -t $cpu -k 30 -R '\@RG\\tID:$name\\tLB:$name\\tSM:$name\\tPL:ILLUMINA' $ref $fq1 $fq2 | awk '\$6 !~ /H/' | samtools view -Sb -t $ref.fai -o $outdir/$name.dh.bam - && echo bwa $name job done\n";
		print BAMLST "$name\t$outdir/$name.dh.bam\n";
	}
	print OUT "echo ends\ndate\n";
	close OUT;
	close BAMLST;

	`chmod +x $outdir/$sh.$name.sh`;
	`sh $outdir/$sh.$name.sh &> $outdir/$sh.$name.sh.log &`;
}
close LIST;
