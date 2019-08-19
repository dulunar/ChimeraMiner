# ChimeraMiner
A pipeline for searching the chimeric sequences in the NGS sequence data.

## Test
The Test folder contains an example. It turns out that all the scripts are running. You can check out how to use them.

## Generate_Shell_Finder.pl	
Generate a shelle script for run the ChimeraMiner.

## Insertion.SRExtract.ReConFastq.pl
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

## SearchOverlapSEchimera.pl	
Search the overlap of the mapped two segment. This only for single-end chimeric sequences.

## A1_Extract.Chimeras.pl
## A2_TransFormat.pl 
## A3_GetInfo.Chimeras.PRE.pl
These 3 scripts provide filtering, classification, and statistics for single-end chimeric sequences.

## the usage of bwa
see the details in the https://github.com/lh3/bwa.

## Contact:
We will be pleased to address any question or concern you may have with the ChimeraMiner: nlu@seu.edu.cn

## Citing ChimeraMiner
If you use ChimeraMiner in your work, please cite:
> Lu, N.; Li, J.; Bi, C.; Guo, J.; Tao, Y.; Luan, K.; Tu, J.; Lu, Z. ChimeraMiner: An Improved Chimeric Read Detection Pipeline and Its Application in Single Cell Sequencing. Int. J. Mol. Sci. 2019, 20, 1953.
