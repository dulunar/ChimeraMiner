#!/usr/bin/perl -w 
use strict;
use Getopt::Long;
my $usage = "
Usage: perl $0 
	-d <string> <the directory of all chimeras file>
	-n <string> <the name of this sample>
	-L <string> <min length of segment> <Default:30>
Example: perl $0 -d /mainsd/luna/work/haplotype_phase/new_pipeline/alignment/SRX247249/chimeras -n SRX247249 -L 20\n";

my ($in,$dir,$name,$minLen,$help);

GetOptions(
	'd=s'   => \$dir,
	'help|?' => \$help,
	'L=s'	=>	\$minLen,
	'n=s'	=>	\$name,
);

#$in ||= "YH_Bulk.ALL.CHR";

die "$usage\n" if ($help || !$dir);

$minLen ||= 30;
`mkdir -p $dir/chimeras_analysis` if(!(-d "$dir/chimeras_analysis"));
my $out = "$dir/chimeras_analysis/normal.txt";
my $out_1 = "$dir/chimeras_analysis/error.txt";

open OUT, ">$out";
open OUT_1, ">$out_1";

for my $chr (1..22,"X","Y","MT"){
	if(-e "$dir/$name.chr$chr"){$in = "$dir/$name.chr$chr";}
	elsif(-e "$dir/$name.chr$chr.gz"){$in = "$dir/$name.chr$chr.gz";}
	else{print "no file in this directory $dir\n";}
	
	next unless $in;

	if($in =~ /\.gz$/){
		open IN,"gzip -dc $in |" || die $!;	
	}
	else{open IN, "<$in" || die $!;}

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
		my ($distance) = $temp_3[5] =~ /(.*)nt/;

		if ($length_2 < $minLen || $length_1 < $minLen) {
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
}
close OUT;
close OUT_1;
