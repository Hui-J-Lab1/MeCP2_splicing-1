SAM1=WT
SAM2=MeCP2_KO
SAM3=RBFOX2_KI
SAM1_1=${SAM1}_1
SAM1_2=${SAM1}_2
SAM2_1=${SAM2}_1
SAM2_2=${SAM2}_2
SAM3_1=${SAM3}_1
SAM3_2=${SAM3}_2

BAM_DIR=/cluster/home/xfu/Project/HJY/20181121JY/RNA-Seq/bam
GTF=/cluster/home/xfu/Gmatic6/gene/hg38_v26/gencode_hg38_v26.gtf

## bam list
echo "$BAM_DIR/$SAM1_1.bam,$BAM_DIR/$SAM1_2.bam" > b1.txt
echo "$BAM_DIR/$SAM2_1.bam,$BAM_DIR/$SAM2_2.bam" > b2.txt
echo "$BAM_DIR/$SAM3_1.bam,$BAM_DIR/$SAM3_2.bam" > b3.txt

mkdir rMATS_out
rmats.py --b1 b1.txt --b2 b2.txt --gtf $GTF --od rMATS_out/${SAM1}_vs_${SAM2} -t paired --nthread 20 --readLength 101 --tstat 20 --libType fr-firststrand 
rmats.py --b1 b1.txt --b2 b3.txt --gtf $GTF --od rMATS_out/${SAM1}_vs_${SAM3} -t paired --nthread 20 --readLength 101 --tstat 20 --libType fr-firststrand 

## select some columns
mkdir rMATS_out_reformat
./script/reformat_rMATS_out.sh ${SAM1}_vs_${SAM2}
./script/reformat_rMATS_out.sh ${SAM1}_vs_${SAM3}

## merge tables
mkdir tables
Rscript script/merge_AS_out.R ${SAM1}_vs_${SAM2} ${SAM1}_vs_${SAM3}
