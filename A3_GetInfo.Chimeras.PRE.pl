#!/usr/bin/perl -w 
use strict;
use Getopt::Long;
my $usage = "

Usage: perl $0 
	-i <string> <in.input filter chimeras file> 
	-d <string> <the valid direct chimeras output file> 
	-v <string> <the valid inverted chimeras output file>
	-L <INT> <the minimum length of segment> <default: 30>
	-s <INT> <the smallest distance of two segment> <default: 25>
	-b <INT> <the biggest distance of two segment> <default: 5000>
\n";

my ($in,$dire,$vert,$len,$smd,$bgd,$help);

$len ||= 30;
$smd ||= 25;
$bgd ||= 5000;


GetOptions(
	'i=s'   => \$in,
	'd=s'   => \$dire,
	'v=s'	=> \$vert,
	'L=i'	=> \$len,
	's=i'	=> \$smd,
	'b=i'	=> \$bgd,
	'help|?' => \$help,
);

die "$usage\n" if ($help || !$in || !$vert || !$dire);

open IN, "<$in" || die $!;
open OD, ">$dire" || die $!;
open OV, ">$vert" || die $1;

my $tag;
my ($invert,$direct) = (0)x2;
my ($La,$Lb,$Lc,$Ld) = (0)x4;
while (my $line_1 = <IN>) {
	my $line_2 = <IN>;
	my $line_3 = <IN>;
	my $line_4 = <IN>;
	my $line_5 = <IN>;
	my $line_6 = <IN>;
	my ($str_1,$seq1) = (split /\s+/, $line_2)[1,4];
	my ($str_2,$seq2) = (split /\s+/, $line_3)[1,4];
	next unless( length($seq1) >= $len );
	next unless( length($seq2) >= $len );
	my ($lap, $dis) = (split /\t/,$line_4)[1,2];
	my $lenp = (split /\s+/,$lap)[0];
	my $distance = (split /\s+/,$dis)[1];
	next unless($lenp ne "" && $lenp >= 3);
	next unless(abs($distance) <= $bgd && abs($distance) >= $smd);
	
	if ($str_1 eq "+" && $str_2 eq "-"){
		$La ++;$tag = "invert";$invert++;
		print OV "$line_1$line_2$line_3$line_4$line_5$line_6";
	}
	elsif ($str_1 eq "-" && $str_2 eq "+"){
		$Lb ++;$tag = "invert";$invert++;
		print OV "$line_1$line_2$line_3$line_4$line_5$line_6";
	}
	elsif ($str_1 eq "+" && $str_2 eq "+"){
		$Lc ++;$tag = "direct";$direct++;
		print OD "$line_1$line_2$line_3$line_4$line_5$line_6";
	}
	elsif ($str_1 eq "-" && $str_2 eq "-"){
		$Ld ++;$tag = "direct";$direct++;
		print OD "$line_1$line_2$line_3$line_4$line_5$line_6";
	}
}

print  "+-\t-+\t++\t--\tInvertChimera\tDirectChimera\n";
print "$La\t$Lb\t$Lc\t$Ld\t$invert\t$direct\n";

close IN;close OV;close OD;
