#!/usr/bin/perl -w 
use strict;
use Getopt::Long;
my $usage = "
Usage: perl $0 
	-i <string> <the chromosome 1 chimeras file>
	-n <string> <the name of this sample>
	-L <string> <min length of segment> <Default:30>
	-c <string> <the valid chimeras output>
	-e <string> <the error chimeras output>

Example: perl $0 -i /home/luna/work/Chimeras/underwrite/sc_mda/statistics/rawchimeras/HUMDA.chr1 -n HUMDA -L 30 -c HUMDA.valid.chimeras -e HUMDA.err
\n";

my ($in,$ochi,$oerr,$name,$minLen,$help);

GetOptions(
	'i=s'	=>	\$in,
	'h|?'	=>	\$help,
	'L=s'	=>	\$minLen,
	'n=s'	=>	\$name,
	'c=s'	=>	\$ochi,
	'e=s'	=>	\$oerr,
);

die "$usage\n" if ($help || !$in);

$minLen ||= 30;

open OUT, ">$ochi";
open OUT_1, ">$oerr";

if($in =~ /\.gz$/){open IN,"gzip -dc $in |" || die $!;}
else{open IN, "$in" || die $!;}

while (my $line_1 = <IN>) {
	my $line_2 = <IN>;
	my $line_3 = <IN>;
	my $line_4 = <IN>;
	my $line_5 = <IN>;
	my $line_6 = <IN>;
	chomp $line_1;
	chomp $line_2;
	chomp $line_3;
	chomp $line_4;
	chomp $line_5;
	chomp $line_6;

	my @temp_1 = split /\t/, $line_3;
	my $str_1 = $temp_1[1];
	my $length_1 = length($temp_1[3]);

	my @temp_2 = split /\t/, $line_4;
	my $str_2 = $temp_2[1];
	my $length_2 = length($temp_2[3]);

	my @temp_3 = split /\s/, $line_5;
	my ($overlap) = $temp_3[3] =~ /(.*)nt/;
	#my ($distance) = $temp_3[5] =~ /(.*)nt/;
	my $distance = $temp_2[2] - $temp_1[4];
	$line_5 = "$temp_3[0]\t$temp_3[1]\t$temp_3[2]\t$temp_3[3]\t$temp_3[4]\t${distance}nt";

	if ($length_2 < $minLen || $length_1 < $minLen || (abs($distance) > 5000) ) {
		print OUT_1 "$line_1\n$line_2\n$line_3\n$line_4\n$line_5\n$line_6\n";
	}
	elsif (($str_1 eq $str_2) && ((abs($distance) <= 2))) {
		print OUT_1 "$line_1\n$line_2\n$line_3\n$line_4\n$line_5\n$line_6\n";
	}
	elsif(abs($overlap) <= 2) {
		print OUT_1 "$line_1\n$line_2\n$line_3\n$line_4\n$line_5\n$line_6\n";
	}
	else {
		print OUT "$line_1\n$line_2\n$line_3\n$line_4\n$line_5\n$line_6\n";
	}
}

close IN;
close OUT;
close OUT_1;
