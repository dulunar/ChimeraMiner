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

## Test
The Test folder contains an example. It turns out that all the scripts are running. You can check out how to use them.

## Generate_Shell_Finder.pl	
Generate a shelle script for run the ChimeraMiner.

## Insertion.SRExtract.ReConFastq.pl
```shell
for first alignment file:
	Step 1:
			seleted first-type chimeric reads	(++ or --).
	Step 2:
			Extract the sequences which are soft-aligned with xxSxxM , xxMxxS or xxSxxMxxS in alignment format from BWA-aligned SAM file.
			        Any Soft-aligned fragment mustn't shorter than log4(Chromosome_Length) nt.
	Step 3:
			Reconstruct FASTQ file from all qualified soft-aligned reads. 
			Qualified reads would be cut into two or three part according to the soft-aligned format.
			All reads would be separated into different chromosomes according to the primary alignment results by BWA.
			Cut reads generated in here would be aligned to Hg19 reference by using BWA-sampe mode.
```

## SearchOverlapSEchimera.pl	
Search the overlap of the mapped two segment. This only for single-end chimeric sequences.

## A1_Extract.Chimeras.pl A2_TransFormat.pl A3_GetInfo.Chimeras.PRE.pl
These 3 scripts provide filtering, classification, and statistics for single-end chimeric sequences.

## the usage of bwa
see the details in the [page](https://github.com/lh3/bwa).

## Contact:
We will be pleased to address any question or concern you may have with the ChimeraMiner: nlu@seu.edu.cn

## Citing ChimeraMiner
If you use ChimeraMiner in your work, please cite:
> Lu, N.; Li, J.; Bi, C.; Guo, J.; Tao, Y.; Luan, K.; Tu, J.; Lu, Z. ChimeraMiner: An Improved Chimeric Read Detection Pipeline and Its Application in Single Cell Sequencing. Int. J. Mol. Sci. 2019, 20, 1953.
