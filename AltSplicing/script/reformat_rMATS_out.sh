VS=$1

cat rMATS_out/$VS/SE.MATS.JCEC.txt   | awk '{print $4":"$5":"$6":"$7":"$8":"$9":"$10":"$11"\t"$2"\t"$3"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}'|sed 's/"//g'|sed 's/chrChr/Chr/' > rMATS_out_reformat/${VS}.SE
cat rMATS_out/$VS/RI.MATS.JCEC.txt   | awk '{print $4":"$5":"$6":"$7":"$8":"$9":"$10":"$11"\t"$2"\t"$3"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}'|sed 's/"//g'|sed 's/chrChr/Chr/' > rMATS_out_reformat/${VS}.RI
cat rMATS_out/$VS/A3SS.MATS.JCEC.txt | awk '{print $4":"$5":"$6":"$7":"$8":"$9":"$10":"$11"\t"$2"\t"$3"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}'|sed 's/"//g'|sed 's/chrChr/Chr/' > rMATS_out_reformat/${VS}.A3SS
cat rMATS_out/$VS/A5SS.MATS.JCEC.txt | awk '{print $4":"$5":"$6":"$7":"$8":"$9":"$10":"$11"\t"$2"\t"$3"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}'|sed 's/"//g'|sed 's/chrChr/Chr/' > rMATS_out_reformat/${VS}.A5SS
cat rMATS_out/$VS/MXE.MATS.JCEC.txt  | awk '{print $4":"$5":"$6":"$7":"$8":"$9":"$10":"$11":"$12":"$13"\t"$2"\t"$3"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25}'|sed 's/"//g'|sed 's/chrChr/Chr/' > rMATS_out_reformat/${VS}.MXE
