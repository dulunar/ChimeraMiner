#########################################################################
# File Name: format.bwa2soap.pl
# Author: Luna
# Mail: nlu@seu.edu.cn
# Created Time: Mon 17 Dec 2018 07:11:16 PM CST
#########################################################################
#!/usr/bin/perl -w
 use Getopt::Long;
 my ($in,$out,$help);

 GetOptions
 (
  "i|in=s"=>\$in,
  "o|out=s"=>\$out,
  "h|help|?"=>\$help,
 );

my $usage=<<INFO;

Usage:

perl $0 [options]
Options:

	-i <file> <in.input file>
	-o <string> <out.output file>

INFO

die $usage if ($help || !$in || !$out);

open IN,"$in" || die $!;
open OUT,"> $out" || die $!;
while(my $a=<IN>){
	my $b=<IN>;
	my $c=<IN>;
	my $d=<IN>;
	my $e=<IN>;
	my $f=<IN>;
	chomp $a;chomp $b;chomp $c;chomp $d;chomp $e;chomp $f;
	my ($chr,$str1,$bg1,$seq1,$ed1) = (split /\s+/,$c)[0,1,2,3,4];
	my ($lhr,$str2,$bg2,$seq2,$ed2) = (split /\s+/,$d)[0,1,2,3,4];
	my ($seq,$lap,$tan) = (split /\s+/,$e)[2,3,5];
	my ($len) = $lap =~ /(.*)nt/;
	my ($dis) = $tan =~ /(.*)nt/;
	my ($read,$id)=(split /\s+/,$f)[0,1];
	print OUT "$b\t$id\n";
	print OUT "$chr\t$str1\t$bg1\t -> $seq1 -> $ed1\n";
	print OUT "$lhr\t$str2\t$bg2\t -> $seq2 -> $ed2\n";
	print OUT "overlap: $seq\t$len nt\tdistance: $dis nt\n";
	print OUT "$read\n\n";
}

close OUT;
close IN;
