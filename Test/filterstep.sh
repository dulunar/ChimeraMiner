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
perl $dir/A1_Extract.Chimeras.pl -d $dir/Test/chimeras -n Test -L 20
perl $dir/A2_TransFormat.pl -i $dir/Test/chimeras/chimeras_analysis/normal.txt -o $dir/Test/chimeras/chimeras_analysis/normal.TransFormat
perl $dir/A3_GetInfo.Chimeras.PRE.pl -i $dir/Test/chimeras/chimeras_analysis/normal.TransFormat -d $dir/Test/chimeras/chimeras_analysis/Test.direct -v $dir/Test/chimeras/chimeras_analysis/Test.inverted  > Test.stat
