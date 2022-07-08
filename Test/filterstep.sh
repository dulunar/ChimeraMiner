#!/bin/bash
#########################################################################
# FileName: filterstep.sh
# Version: 23178b66-4cec-4e30-a952-b4f77fdb7c8b
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Tue 14 Apr 2020 05:17:28 PM CST
#########################################################################

dir=/home/luna/Desktop/Tools/git/ChimeraMiner

cd $dir/Test
# delete the log file generated from workstep.sh
rm -rf run.4Searcha* log.step4.chr*

# Filtering, classification, and statistics, All the parameters represent what can be seen in the script.
perl ChimerasDownstream.pl -d $dir/Test/chimeras -n Test -L 30 -s 25 -b 5000 > $dir/Test/Test.chimera.stat.txt 
