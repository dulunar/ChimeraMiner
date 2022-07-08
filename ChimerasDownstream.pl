#!/usr/bin/perl -w 
use strict;

my $red = "\033[0;31m";
my $rend = "\033[0m";

use Getopt::Long;
my $usage = "
Usage: perl $0
   For Extracting, transformating, and getinfo of chimeras.
	-d <string> <the directory of all chimeras file>
	-n <string> <the name of this sample>
	-L <string> <min length of segment> <Default:30>
	-s <INT> <the smallest distance of two segment> <default: 25>
	-b <INT> <the biggest distance of two segment> <default: 5000>

   Example: 
   ${red}perl $0 -d /mainsd/luna/work/haplotype_phase/new_pipeline/alignment/SRX247249/chimeras -n SRX247249 -L 30 -s 25 -b 5000
\n";

my ($in, $dir, $name, $minLen, $help, $smd, $bgd);

GetOptions(
	'd=s'		=>	\$dir,
	'help|?'	=>	\$help,
	'L=s'		=>	\$minLen,
	'n=s'		=>	\$name,
	's=i'		=>	\$smd,
	'b=i'		=>	\$bgd,
);
die "$usage\n" if ($help || !$dir || !$name);

$minLen ||= 30;
$smd ||= 25;
$bgd ||= 5000;

`mkdir -p $dir/chimeras_analysis` if(!(-d "$dir/chimeras_analysis"));
my $out = "$dir/chimeras_analysis/normal.txt";
my $out_1 = "$dir/chimeras_analysis/error.txt";

open OUT, ">$out";
open OUT_1, ">$out_1";

open OT,"> $dir/chimeras_analysis/$name.TransFormat" || die $!;
open OD,"> $dir/chimeras_analysis/$name.direct" || die $!;
open OV,"> $dir/chimeras_analysis/$name.inverted" || die $!;

my ($invert,$direct) = (0)x2;
my ($La,$Lb,$Lc,$Ld) = (0)x4;
my $tag;

die $red."$dir not exists $rend\n" if(! -e $dir);

for my $chr (1..22,"X","Y", "MT") {
	if(-e "$dir/$name.chr$chr") {
		$in = "$dir/$name.chr$chr";
	}
	elsif(-e "$dir/$name.chr$chr.gz") {
		$in = "$dir/$name.chr$chr.gz";
	}
	else{
		die "${red}no $chr file in this directory $dir$rend\n";
	}
	
	next unless $in;

	if($in =~ /\.gz$/){
		open IN,"gzip -dc $in |" || die $!;	
	}
	else{
		open IN, "$in" || die $!;
	}

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

		my ($chr,$str1,$bg1,$seq1,$ed1) = (split /\t/,$line_3)[0,1,2,3,4];
		my ($lhr,$str2,$bg2,$seq2,$ed2) = (split /\t/,$line_4)[0,1,2,3,4];
		my ($seqp, $lap, $tan) = (split /\t/,$line_5)[2,3,5];
		my ($lenp) = $lap =~ /(.*)nt/;
		my ($dis) = $tan =~ /(.*)nt/;
		my ($read, $id)=(split /\t/,$line_6)[0,1];
		
		my $trans = "$line_2\t$id\n$chr\t$str1\t$bg1\t -> $seq1 -> $ed1\n$lhr\t$str2\t$bg2\t -> $seq2 -> $ed2\noverlap: $seqp\t$lenp nt\tdistance: $dis nt\n$read\n\n";

		my $len1 = length($seq1);

		my $len2 = length($seq2);
		
		if ($len2 < $minLen || $len1 < $minLen) {
			print OUT_1 "$line_1\n$line_2\n$line_3\n$line_4\n$line_5\n$line_6\n";
		}
		elsif(abs($dis) > $bgd && abs($dis) > $smd) {
			print OUT_1 "$line_1\n$line_2\n$line_3\n$line_4\n$line_5\n$line_6\n";
		}
		elsif(abs($lenp) <= 2 || $lenp eq "") {
			print OUT_1 "$line_1\n$line_2\n$line_3\n$line_4\n$line_5\n$line_6\n";
		}
		else{
			print OUT "$line_1\n$line_2\n$line_3\n$line_4\n$line_5\n$line_6\n";
			print OT "$trans";
		
			if ($str1 eq "+" && $str2 eq "-"){
				$La ++;$tag = "invert";$invert++;
				print OV "$trans";
			}
			elsif ($str1 eq "-" && $str2 eq "+"){
				$Lb ++;$tag = "invert";$invert++;
				print OV "$trans"; 
			}
			elsif ($str1 eq "+" && $str2 eq "+"){
				$Lc ++;$tag = "direct";$direct++;
				print OD "$trans";
			}
			elsif ($str1 eq "-" && $str2 eq "-"){
				$Ld ++;$tag = "direct";$direct++;
				print OD "$trans";
			}
		}

	}
	close IN;
}
close OUT;
close OUT_1;
close OT;
close OV;
close OD;

print "Title\t+-\t-+\t++\t--\tInvertChimera\tDirectChimera\n";
print "$name\t$La\t$Lb\t$Lc\t$Ld\t$invert\t$direct\n";
