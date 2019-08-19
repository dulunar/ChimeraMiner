echo starts
date
perl /home/luna/work/Chimeras/ChimeraMinerInsertion.SRExtract.ReConFastq.pl -i /home/luna/work/Chimeras/ChimeraMiner/Test/Test.dh.bam -m Test -d /home/luna/work/Chimeras/ChimeraMiner/Test && echo first split and find INSERT chimerasjob done
date
sh /home/luna/work/Chimeras/ChimeraMiner/Test/Chr_split/run.aln.sh &> /home/luna/work/Chimeras/ChimeraMiner/Test/Chr_split/run.aln.sh.log && echo realign soft-alignment reads to reference
split -l 6 /home/luna/work/Chimeras/ChimeraMiner/Test/run.Test.4Search.sh /home/luna/work/Chimeras/ChimeraMiner/Test/run.4Search
chmod +x /home/luna/work/Chimeras/ChimeraMiner/Test/run.4Search*
ls /home/luna/work/Chimeras/ChimeraMiner/Test/run.4Search* | perl -ne 'chomp;`nohup sh $_ &> $_.log &`;'