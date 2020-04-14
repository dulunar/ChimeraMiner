#!/bin/bash
#########################################################################
# FileName: workstep.sh
# Version: 7c6df8b9-7e6c-43c5-8fd4-75c1e7a5f0b8
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Sun 18 Aug 2019 11:00:56 PM CST
#########################################################################

# PS: All results were put in this folder.
dir=/home/luna/Desktop/Tools/git/ChimeraMiner
if [ ! -f "$dir/Test" ]; then
	mkdir -p $dir/Test 
fi

cd $dir/Test 

ref=/home/luna/Desktop/database/homo_bwa/hsa.fa

# step 1, run bwa mem, align the reads to the reference， here we used the hg19 which downloaded from UCSC.
bwa mem -t 20 -k 30 -R '@RG\tID:Test\tLB:Test\tSM:Test\tPL:ILLUMINA' $ref test_R1.fq.gz test_R2.fq.gz | awk '$6 !~ /H/' | samtools view -Sb -t ${ref}.fai -o $dir/Test/Test.dh.bam - && echo bwa Test job done

echo -e "Test\t$dir/Test/Test.dh.bam" > bam.lst

# step 2, generate a shell script for run ChimeraMiner, in here "-L Options" means the min Length of each segment.
perl $dir/Generate_Shell_Finder.pl -i bam.lst -o runFinder.Test.sh -L 20 -r $ref

# step3, run the shell script, when search the overlap, We carry out each chromosome independently.
sh runFinder.Test.sh
