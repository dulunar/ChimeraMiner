# ChimeraMiner
A pipeline for searching the chimeric sequences in the NGS sequence data.

## Reference Genome Preparation
In my server, my reference genome fold is:

/home/luna/Desktop/database/homo_bwa

In this fold, we have 26 fasta files (chr{1..22, X, Y, MT}.fa hsa.fa), each fasta file have indexed with bwa index.

In additional, `hsa.fa` include all sequence from all chromosome (chr{1..22, X, Y, MT}).

Index for each fasta file:

```shell
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
do
	bwa index -a bwtsw chr$chr.fa &>> bwa-index.log
	samtools faidx chr$chr.fa &>> samtool-index.log
done
bwa index -a bwtsw hsa.fa &>> bwa-index.log
samtools faidx hsa.fa &>> samtool-index.log
```
## Perl Modules Preparation
In this pipeline, we need some Perl Modules, we should install these modules before run it.
1. Getopt::Long
2. Cwd qw(abs_path)
3. File::Basename
4. File::Spec

I recommand to use "cpanm" to install Perl Modules.

## Test
[The Test folder contains an example](https://github.com/dulunar/ChimeraMiner/tree/master/Test). It turns out that all the scripts are running. You can check out how to use the pipeline. 
In this folder, just run workstep.sh first, this shell will generate bam file and chimera's files. When all works in workstep.sh finished, then run filterstep.sh, this shell will deal with the chimera's files and count.

### First, align reads to reference
use bwa mem to map the reads to hg19, generate a list contains "SampleID	BAMFile", each sample each line;
```shell
ref=/home/luna/Desktop/database/homo_bwa/hsa.fa
bwa mem -t 20 -k 30 -R '@RG\tID:Test\tLB:Test\tSM:Test\tPL:ILLUMINA' $ref test_R1.fq.gz test_R2.fq.gz | awk '$6 !~ /H/' | samtools view -Sb -t ${ref}.fai -o $dir/Test/Test.dh.bam -
echo -e "Test\t$dir/Test/Test.dh.bam" > bam.lst
```

###  Generate_Shell_Finder.pl	
Use bam.lst as input file. Generate a shelle script for running the ChimeraMiner.
```shell
perl Generate_Shell_Finder.pl -i bam.lst -o runFinder.Test.sh -L 20 -r $ref
```

### runFinder.Test.sh
run the shell script, do some works:
1. Insertion.SRExtract.ReConFastq.pl for searching insertion chimeras, extracting soft-clipped alignment reads as candidate single-ended chimeric reads (direct and inverted) and re-constructring  pe fastqs (each chromosome has two fastq files) for candidate chimeras
2. aligned chromosome fastq files to chromosome reference respectively.
3. SearchOverlapSEchimera.pl for searching the overlap sequence of two adjacent segments of  candidate single-ended chimeric reads, when searching the overlap sequence, we carry out each chromosome independently.

### ChimerasDownstream.pl
This script is used for downstream analysis of chimeras:
1. Merge chimeras of each chromosome to a file
2. Tranform the raw format of chimeras to a better format for viewing
3. Extracting the direct chimera and inverted chimera to different files.
4. Count the chimera types'++' '--' '+-' '-+' , and the number of direct chimera and inverted chimera.

## the usage of bwa
see the details in the [page](https://github.com/lh3/bwa).

## Contact:
We will be pleased to address any question or concern you may have with the ChimeraMiner: nlu@seu.edu.cn

## Citing ChimeraMiner
If you use ChimeraMiner in your work, please cite:
> Lu, N.; Li, J.; Bi, C.; Guo, J.; Tao, Y.; Luan, K.; Tu, J.; Lu, Z. ChimeraMiner: An Improved Chimeric Read Detection Pipeline and Its Application in Single Cell Sequencing. Int. J. Mol. Sci. 2019, 20, 1953.
