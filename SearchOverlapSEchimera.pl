#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $usage = "
for first alignment file:
1.seleted first-type chimeric reads	(++ or --).
2.search the overlap of the mapped two segment.
Usage: perl $0 
	-i <bam file>
	-o <out file>
	-d <output directory>
	-L <min Length of each segment> <default:30>
";

my ($in,$out,$dir,$minLen,$help);

GetOptions(
		'i=s'	=>	\$in,
		'o=s'	=>	\$out,
		'd=s'	=>	\$dir,
		"L=s"	=>	\$minLen,
		'h|?'	=>	\$help,
		);

die "$usage\n" if ($help || !$in || !$out || !$dir);

$minLen ||= 30;

`mkdir -p $dir/chimeras`if(!(-d "$dir/chimeras"));
my $samtools = "/home/luna/Desktop/Software/samtools/samtools";

if($in =~ /\.bam$/){
	open IN, "$samtools view $in | awk '\$6 !~ /H/' |" or die "$!";
}
elsif($in =~ /\.gz$/){
	open IN, "gzip -dc $in | awk '\$6 !~ /H/' |" or die "$!";
}
else{
	open IN, "cat $in | awk '\$6 !~ /H/' |" or die "$!";
}

open OUT, "> $dir/chimeras/$out" || die $!;
open ERR, "> $dir/chimeras/$out.err" || die $!;
while(<IN>){
	chomp;
	next unless($_ !~ /^@/);
	my $line = <IN>;
	chomp $line;
	my ($id_1, $flag_1, $chr_1, $pos_1, $qua_map_1,$format_1, $seq_1, $qua_1) = (split /\t/, $_)[0, 1, 2, 3, 4, 5, 9, 10];
	my ($id_2, $flag_2, $chr_2, $pos_2, $qua_map_2, $format_2, $seq_2, $qua_2) = (split /\t/, $line)[0, 1, 2, 3, 4, 5, 9, 10];

	my $str_1 = strand_code($flag_1);
	my $str_2 = strand_code($flag_2);
	chomp $seq_1;chomp $seq_2;	
	$seq_1 = seq_trans($seq_1) if($str_1 eq "-");
	$seq_2 = seq_trans($seq_2) if($str_2 eq "-");
	my ($id_real) = $id_1 =~ /(.*)\.\w+$/;
	my $map1 = length $seq_1;
	my $map2 = length $seq_2;

	if($format_1=~/S/ || $format_2=~/S/ || $format_1 eq "*" || $format_2 eq "*" || $map1 < $minLen || $map2 < $minLen){
		print ERR "$_\n";
		print ERR "$line\n";
	}
	else{
		my ($I_1, $D_1, $I_2, $D_2) = (0, 0, 0, 0);
		my $pos_1_start = $pos_1;	
		($I_1) = $format_1 =~ /(\d+)I/ if($format_1 =~ /I/);		($D_1) = $format_1 =~ /(\d+)D/ if($format_1 =~ /D/);
		my $pos_1_end = $pos_1 + $map1 + $D_1 - $I_1 -1;

		my $pos_2_start = $pos_2;	
		($I_2) = $format_2 =~ /(\d+)I/ if($format_2 =~ /I/);		($D_2) = $format_2 =~ /(\d+)D/ if($format_2 =~ /D/);
		my $pos_2_end = $pos_2 + $map2 + $D_1 - $I_1 - 1;

		my $segment_ref_sequence;	my $reverse_extend_forward;		my $reverse_extend_end;
		my $tag = "D"; #解释overlap处于reads的哪一segment上,F代表前一段，E代表后一段；

			my $cut_forward_seq;	#切割segment1的后端;
		my $cut_end_seq;	#切割segment2的前端;

		my $length_reads = $map1 + $map2;
		my $ref_reads_sequence;
		my ($overlap_forward, $overlap_forward_length, $overlap_end, $overlap_end_length, $overlap_seq, $overlap_seq_length) = ("", 0, "", 0, "", 0);

		if($id_1 =~ /\.M$/ || $id_1 =~ /\.2$/ || $id_1 =~ /\.S$/ || $id_1 =~ /\.1$/){
			my $tmp_reads_referecne_end = $length_reads + $pos_1_start - 1;
			$ref_reads_sequence = ref_reads_seq($chr_1, $pos_1_start, $tmp_reads_referecne_end);
			$cut_forward_seq = substr($seq_1, -31);   ##切割segment1的后端31bp;
			$cut_forward_seq = length_cut_1($cut_forward_seq);
			$cut_end_seq = substr($seq_2, 0, 31);     ##切割segment2的前端31bp;
			$cut_end_seq = length_cut_2($cut_end_seq);

			if($str_1 eq "+"){
				my $tmp_start = $pos_1_end + 1;
				my $tmp_end = $pos_1_end + 31;
				$segment_ref_sequence = ref_reads_seq($chr_1, $tmp_start, $tmp_end);
				$reverse_extend_forward = length_cut_2(uc($segment_ref_sequence));
				($overlap_end, $overlap_end_length) = compare_sequence_end($cut_end_seq, $reverse_extend_forward);

				if($str_2 eq "+"){
					$tmp_start = $pos_2_start - 31;
					$tmp_end = $pos_2_start - 1;
					$segment_ref_sequence = ref_reads_seq($chr_2, $tmp_start, $tmp_end);
					$reverse_extend_end = length_cut_1(uc($segment_ref_sequence));
					($overlap_forward, $overlap_forward_length) = compare_sequence_forward($cut_forward_seq, $reverse_extend_end);
					($overlap_forward_length >= $overlap_end_length)?(($overlap_seq = $overlap_forward) && ($tag = "F")):(($overlap_seq = $overlap_end) && ($tag = "E"));
					$overlap_seq_length = length $overlap_seq;
					my $distance = $pos_2_start - $pos_1_end;
					print OUT "##hg19\t$chr_1\t+\t$pos_1_start\t$tmp_reads_referecne_end\n$ref_reads_sequence\n";
					print OUT "$chr_1\t$str_1\t$pos_1_start\t$seq_1\t$pos_1_end\n";
					print OUT "\t$str_2\t$pos_2_start\t$seq_2\t$pos_2_end\n";
					print OUT "overlap:\t$tag\t$overlap_seq\t${overlap_seq_length}nt\tdistance:\t${distance}nt\n";
					print OUT "$seq_1$seq_2\t$id_real\n";
				}
				else{
					$tmp_end = $pos_2_end + 31;
					$tmp_start = $pos_2_end + 1;
					$segment_ref_sequence = ref_reads_seq($chr_2, $tmp_start, $tmp_end);
					$reverse_extend_end = length_cut_1(seq_trans(uc($segment_ref_sequence)));
					($overlap_forward, $overlap_forward_length) = compare_sequence_forward($cut_forward_seq, $reverse_extend_end);
					($overlap_forward_length >= $overlap_end_length)?(($overlap_seq = $overlap_forward) && ($tag = "F")):(($overlap_seq = $overlap_end) && ($tag = "E"));
					$overlap_seq_length = length $overlap_seq;
					my $distance = $pos_2_start - $pos_1_end - $overlap_seq_length;
					print OUT "##hg19\t$chr_1\t+\t$pos_1_start\t$tmp_reads_referecne_end\n$ref_reads_sequence\n";
					print OUT "$chr_1\t$str_1\t$pos_1_start\t$seq_1\t$pos_1_end\n";
					print OUT "\t$str_2\t$pos_2_end\t$seq_2\t$pos_2_start\n";
					print OUT "overlap:\t$tag\t$overlap_seq\t${overlap_seq_length}nt\tdistance:\t${distance}nt\n";
					print OUT "$seq_1$seq_2\t$id_real\n";
				}
			}
			else{
				my $tmp_end = $pos_1_start - 1;
				my $tmp_start = $pos_1_start - 31;
				$segment_ref_sequence = ref_reads_seq($chr_1, $tmp_start, $tmp_end);
				$reverse_extend_forward = length_cut_2(seq_trans(uc($segment_ref_sequence)));
				($overlap_end, $overlap_end_length) = compare_sequence_end($cut_end_seq, $reverse_extend_forward);

				if($str_2 eq "+"){
					$tmp_start = $pos_2_start - 31;
					$tmp_end = $pos_2_start - 1;
					$segment_ref_sequence = ref_reads_seq($chr_2, $tmp_start, $tmp_end);
					$reverse_extend_end = length_cut_1(uc($segment_ref_sequence));
					($overlap_forward, $overlap_forward_length) = compare_sequence_forward($cut_forward_seq, $reverse_extend_end);
					($overlap_forward_length >= $overlap_end_length)?(($overlap_seq = $overlap_forward) && ($tag = "F")):(($overlap_seq = $overlap_end) && ($tag = "E"));
					$overlap_seq_length = length $overlap_seq;
					my $distance = $pos_2_start - $pos_1_end - $overlap_seq_length;
					print OUT "##hg19\t$chr_1\t+\t$pos_1_start\t$tmp_reads_referecne_end\n$ref_reads_sequence\n";
					print OUT "$chr_1\t$str_1\t$pos_1_end\t$seq_1\t$pos_1_start\n";
					print OUT "\t$str_2\t$pos_2_start\t$seq_2\t$pos_2_end\n";
					print OUT "overlap:\t$tag\t$overlap_seq\t${overlap_seq_length}nt\tdistance:\t${distance}nt\n";
					print OUT "$seq_1$seq_2\t$id_real\n";
				}
				else{
					$tmp_end = $pos_2_end + 31;
					$tmp_start = $pos_2_end + 1;
					$segment_ref_sequence = ref_reads_seq($chr_2, $tmp_start, $tmp_end);
					$reverse_extend_end = length_cut_1(seq_trans(uc($segment_ref_sequence)));
					($overlap_forward, $overlap_forward_length) = compare_sequence_forward($cut_forward_seq, $reverse_extend_end);
					($overlap_forward_length >= $overlap_end_length)?(($overlap_seq = $overlap_forward) && ($tag = "F")):(($overlap_seq = $overlap_end) && ($tag = "E"));
					$overlap_seq_length = length $overlap_seq;
					my $distance = $pos_1_start - $pos_2_end;
					print OUT "##hg19\t$chr_1\t+\t$pos_1_start\t$tmp_reads_referecne_end\n$ref_reads_sequence\n";
					print OUT "$chr_1\t$str_1\t$pos_1_end\t$seq_1\t$pos_1_start\n";
					print OUT "\t$str_2\t$pos_2_end\t$seq_2\t$pos_2_start\n";
					print OUT "overlap:\t$tag\t$overlap_seq\t${overlap_seq_length}nt\tdistance:\t${distance}nt\n";
					print OUT "$seq_1$seq_2\t$id_real\n";
				}
			}
		}
		else{print "$_\n$line\n";}
	}
}
close IN;
close OUT;
close ERR;

##根据FLAG值判断reads比对到基因组的正链还是负链;
sub strand_code{
	my $flag = $_[0];#数组是@_, 但是单个参数是$_[0]开始;
	chomp $flag;
	$flag = reverse( sprintf("%b",$flag) );
	my $strand;
	if($flag == 0){
		$strand = "+";
	}
	else{
		my @tmp = split //,$flag;
		$strand = ($tmp[4] == 0)?("+"):("-");
	}
	return $strand;
}

##补充切割得到的sequence的长度到定长;
sub length_cut_1{
	my $seq_cut = $_[0];
	if(length($seq_cut) < 31){
		my $tmp_seq = "N"x(31 - length($seq_cut));
		$seq_cut = $tmp_seq.$seq_cut;
	}
	$seq_cut = uc($seq_cut);
	return $seq_cut;
}

sub length_cut_2{
	my $seq_cut = $_[0];
	if(length($seq_cut) < 31){
		my $tmp_seq = "N"x(31 - length($seq_cut));
		$seq_cut = $seq_cut.$tmp_seq;
	}
	$seq_cut = uc($seq_cut);
	return $seq_cut;
}


##获得全长reads在reference上的序列;
sub ref_reads_seq{
	my ($chr,$start,$end) = @_;
	my $ref_reads_sequence = `$samtools faidx /home/luna/Desktop/database/homo_bwa/$chr.fa $chr:$start-$end`;
	my @ALL = split /\n/,$ref_reads_sequence; 
	$ref_reads_sequence=""; 
	for(my $i = 1;$i<=$#ALL;$i++){
		$ref_reads_sequence .= $ALL[$i];
	}
	$ref_reads_sequence = uc($ref_reads_sequence);
	return $ref_reads_sequence;
}

##对序列进行反向互补转换;
sub seq_trans {
	my $seq = $_[0];
	chomp $seq;
	$seq = uc($seq);
	my $reverse_seq;
	$seq=~tr/ATCG/TAGC/;
	$reverse_seq = reverse $seq;
	return $reverse_seq;
}

##两条序列进行比较得到最大的overlap sequence,允许一个mismatch，但是mismatch不能在overlap的开始位置;
sub compare_sequence_forward{
	my ($cut_forward_seq, $reverse_extend_end) = @_;
	my $overlap = ""; 
	my $length_overlap = 0;
	my @cut_read_seq = split //,$cut_forward_seq;
	my @extend_ref_seq = split //,$reverse_extend_end;
	my $fault = 0;
	my $num = $#cut_read_seq;
	while($num >= 0){
		if($cut_read_seq[$num] ne $extend_ref_seq[$num]){$fault ++;}
		if(($fault == 2) && ($cut_read_seq[$num+1] ne $extend_ref_seq[$num+1])){
			$overlap = substr($overlap, 1);
		}
		last if $fault == 2;
		$overlap = $cut_read_seq[$num].$overlap;
		$num --;
	}
	$length_overlap = length $overlap;
	#print "forward\t$cut_forward_seq\t$reverse_extend_end\t$overlap\t$length_overlap\n";
	return($overlap,$length_overlap);
}

sub compare_sequence_end{
	my ($cut_end_seq, $reverse_extend_forward) = @_;
	my $overlap = "";
	my $length_overlap = 0;
	my @cut_read_seq = split //,$cut_end_seq;
	my @extend_ref_seq = split //,$reverse_extend_forward;
	my $fault = 0;
	my $num = 0;
	while($num <= $#cut_read_seq){
		if($cut_read_seq[$num] ne $extend_ref_seq[$num]){$fault ++;}
		if(($fault == 2) && ($cut_read_seq[$num-1] ne $extend_ref_seq[$num-1])){
			$overlap = substr($overlap, 0, -1);
		}
		last if $fault == 2;
		$overlap .= $cut_read_seq[$num];
		$num ++;
	}
	$length_overlap = length $overlap;
	#print "end\t$cut_end_seq\t$reverse_extend_forward\t$overlap\t$length_overlap\n";
	return($overlap,$length_overlap);

}
