#########################################################################
# FileName: workstep.sh
# Version: 7c6df8b9-7e6c-43c5-8fd4-75c1e7a5f0b8
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Sun 18 Aug 2019 11:00:56 PM CST
#########################################################################

#!/bin/bash

# PS: All results were put in this folder.
mkdir -p /home/luna/work/Chimeras/ChimeraMiner/Test 
cd /home/luna/work/Chimeras/ChimeraMiner/Test  

# step 1, run bwa mem, align the reads to the reference， here we used the hg19 which downloaded from UCSC.
bwa mem -t 20 -k 30 -R '@RG\tID:Test\tLB:Test\tSM:Test\tPL:ILLUMINA' /referenceforbwa/hsa.fa test_R1.fq.gz test_R2.fq.gz | awk '$6 !~ /H/' | samtools view -Sb -t /referenceforbwa/hsa.fa.fai -o /home/luna/work/Chimeras/ChimeraMiner/Test/Test.dh.bam - && echo bwa Test job done

echo -e "Test\t/home/luna/work/Chimeras/ChimeraMiner/Test/Test.dh.bam" > bam.lst

# step 2, generate a shell script for run ChimeraMiner, in here "-L Options" means the min Length of each segment.
perl //home/luna/work/Chimeras/ChimeraMiner/Generate_Shell_Finder.pl -i bam.lst -o runFinder.Test.sh -L 20

# step3, run the shell script, when search the overlap, We carry out each chromosome independently.
sh runFinder.Test.sh

#step 4, Filtering, classification, and statistics, All the parameters represent what can be seen in the script.
perl //home/luna/work/Chimeras/ChimeraMiner/A1_Extract.Chimeras.pl -d /home/luna/work/Chimeras/ChimeraMiner/Test/chimeras -n Test -L 20 
perl //home/luna/work/Chimeras/ChimeraMiner/A2_TransFormat.pl -i normal.txt -o normal.TransFormat
perl //home/luna/work/Chimeras/ChimeraMiner/A3_GetInfo.Chimeras.PRE.pl -i normal.TransFormat -d Test.direct -v Test.inverted  > Test.stat