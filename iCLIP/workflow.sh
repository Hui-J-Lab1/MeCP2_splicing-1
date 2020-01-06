DATA=/var/data/raw_data/HJY/HJY_20190123_YJiang_RBFOX2_iClip_2sam_3rep/Undetermined_S0_L006_R1_001.fastq.gz
barcodeLength=15
minBaseQuality=10
minReadLength=15
maxReadLength=150
adapter=TGAGATCGGAAGAGCGGTTCAG
genomeMappingIndex=/home/xfu/Gmatic7/genome/human/STAR/
GTF=/home/xfu/Gmatic7/gene/human/GRCh38_v29.gtf
FAI=/home/xfu/Gmatic7/genome/human/GRCh38.fa.fai

##====================
## Quality control
##====================

###--------------------------
### General quality check
###--------------------------

#### FastQC on the full data set
mkdir fastqc
fastqc --extract --nogroup --outdir fastqc $DATA

###------------------------------------------
### Quality filter on the barcode regions
###------------------------------------------

#### List of read IDs of reads with high quality barcode regions (using FASTX-Toolkit)
mkdir tmp
zcat $DATA | fastx_trimmer -l $barcodeLength | fastq_quality_filter -q $minBaseQuality -p 100 | awk 'FNR%4==1 { print $1 }' | sed 's/@//' > tmp/data.qualFilteredIDs.list

#### Extract reads of given read IDs (using seqtk) and remove problematic characters and whitespaces from read IDs
./seqtk_batch.sh

###------------------------
### Barcode frequencies 
###------------------------

#### Extract all detected experimental barcodes and their frequencies (x = length of UMI1, y = length of the experimental barcodes)
mkdir results
./barcode_freq_batch.sh

##=================================================
## Demultiplexing, adapter and barcode trimming
##=================================================

#### Demultiplexing, adapter and barcode trimming using Flexbar
mkdir demultiplex
cd demultiplex
../demultiplex_batch.sh

cat x*_Ctrl.fastq.gz > Ctrl.fastq.gz
cat x*_KO_rep1.fastq.gz > KO_rep1.fastq.gz
cat x*_KO_rep2.fastq.gz > KO_rep2.fastq.gz
cat x*_KO_rep3.fastq.gz > KO_rep3.fastq.gz
cat x*_WT_rep1.fastq.gz > WT_rep1.fastq.gz
cat x*_WT_rep2.fastq.gz > WT_rep2.fastq.gz
cat x*_WT_rep3.fastq.gz > WT_rep3.fastq.gz

cd ..
cat demultiplex/*.lengthdist|sort -n|uniq|groupBy -g 1 -c 2 -o sum|sed 's/\./Count/' > results/all_reads.lengthdist

##====================
## Genomic mapping
##====================

#### Genomic mapping using STAR
mkdir mapping
cd mapping
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix Ctrl_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/Ctrl.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix KO_rep1_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/KO_rep1.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix KO_rep2_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/KO_rep2.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix KO_rep3_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/KO_rep3.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix WT_rep1_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/WT_rep1.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix WT_rep2_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/WT_rep2.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix WT_rep3_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/WT_rep3.fastq.gz
find *.bam|parallel --gnu "samtools index {}"
cd ..


##======================================
## Duplicate removal (deduplication)
##======================================

#### Duplicate removal (deduplication) using UMI-tools
mkdir dedup
find mapping/*.bam|sed 's/mapping\///;s/\.bam//'|xargs -I {} umi_tools dedup -I mapping/{}.bam -L dedup/{}.duprm.log -S dedup/{}.duprm.bam --extract-umi-method read_id --method unique


##==========================================
## Extraction of crosslinked nucleotides
##==========================================

#### Convert all read locations to intervals in bed file format using BEDTools
mkdir bed
find dedup/*.bam|sed 's/dedup\///;s/\.bam//'|parallel --gnu "bedtools bamtobed -i dedup/{}.bam > bed/{}.bed"
 
#### Shift intervals depending on the strand by 1 bp upstream using BEDTools
find bed/*.bed|sed 's/bed\///;s/\.bed//'|parallel --gnu "bedtools shift -m 1 -p -1 -i bed/{}.bed -g $FAI > bed/{}.shifted.bed"
 
#### Extract the 5' end of the shifted intervals and pile up into coverage track in bedgraph file format (separately for each strand) using BEDTools (in case of RPM-normalised coverage tracks, use additional parameter -scale with 1,000,000/#mappedReads)
find bed/*.shifted.bed|sed 's/bed\///;s/\.bed//'|parallel --gnu "bedtools genomecov -bg -strand + -5 -i bed/{}.bed -g $FAI > bed/{}.plus.bedgraph"
find bed/*.shifted.bed|sed 's/bed\///;s/\.bed//'|parallel --gnu "bedtools genomecov -bg -strand - -5 -i bed/{}.bed -g $FAI > bed/{}.minus.bedgraph"

#### make tdf file
mkdir track
find bed/*.bedgraph|sed 's/bed\///;s/\.bedgraph//'|parallel --gnu "igvtools toTDF bed/{}.bedgraph bed/{}.tdf /home/xfu/igv/genomes/hg38.genome"


#### Optional convertion of bedgraph files to bw file format files using bedGraphToBigWig of the kentUtils suite
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod +x bedGraphToBigWig 

export LC_COLLATE=C
find bed/*.bedgraph|sed 's/\.bedgraph//'|parallel --gnu "sort -k1,1 -k2,2n {}.bedgraph > {}.sorted.bedgraph"
find bed/*.sorted.bedgraph|sed 's/\.bedgraph//'|parallel --gnu "./bedGraphToBigWig {}.bedgraph $FAI {}.bw"

##========================================================
## crosslink site and deduplication statistics
##========================================================
./Xlink_stat.sh Ctrl
./Xlink_stat.sh KO_rep1
./Xlink_stat.sh KO_rep2
./Xlink_stat.sh KO_rep3
./Xlink_stat.sh WT_rep1
./Xlink_stat.sh WT_rep2
./Xlink_stat.sh WT_rep3

./deduplication_stat.sh Ctrl
./deduplication_stat.sh KO_rep1
./deduplication_stat.sh KO_rep2
./deduplication_stat.sh KO_rep3
./deduplication_stat.sh WT_rep1
./deduplication_stat.sh WT_rep2
./deduplication_stat.sh WT_rep3

# make_stat_table
paste results/dedup_stat*.tsv|cut -f1,2,4,6,8,10,12,14,16 > stat_table.tsv
paste results/Xlink_stat*.tsv|cut -f1,2,4,6,8,10,12,14,16|grep -v 'Ctrl' >> stat_table.tsv

##===============================
## Peak calling with PureCLIP
##===============================

#### Merge BAM files (sampleX.duprm.bam)
mkdir peak
cd peak

samtools merge -f WT_merged.bam ../dedup/WT_rep1_Aligned.sortedByCoord.out.duprm.bam ../dedup/WT_rep2_Aligned.sortedByCoord.out.duprm.bam ../dedup/WT_rep3_Aligned.sortedByCoord.out.duprm.bam
samtools index WT_merged.bam

samtools merge -f KO_merged.bam ../dedup/KO_rep1_Aligned.sortedByCoord.out.duprm.bam ../dedup/KO_rep2_Aligned.sortedByCoord.out.duprm.bam ../dedup/KO_rep3_Aligned.sortedByCoord.out.duprm.bam
samtools index KO_merged.bam

#### Run PureCLIP
GENOME=/home/xfu/Gmatic7/genome/human/GRCh38.fa
pureclip -i KO_merged.bam -bai KO_merged.bam.bai -g $GENOME -ld -nt 10 -o KO_PureCLIP.crosslink_sites.bed -or KO_PureCLIP.crosslink_regions.bed
pureclip -i WT_merged.bam -bai WT_merged.bam.bai -g $GENOME -ld -nt 10 -o WT_PureCLIP.crosslink_sites.bed -or WT_PureCLIP.crosslink_regions.bed

#### Remove 7th column of PureCLIP output file
cat KO_PureCLIP.crosslink_sites.bed | cut -f 1,2,3,4,5,6|grep '^chr' > KO_PureCLIP.crosslink_sites_short.bed
cat WT_PureCLIP.crosslink_sites.bed | cut -f 1,2,3,4,5,6|grep '^chr' > WT_PureCLIP.crosslink_sites_short.bed

bedtools bamtobed -i peak/KO_merged.bam > peak/KO_merged.bed
bedtools bamtobed -i peak/WT_merged.bam > peak/WT_merged.bed

bedtools shift -m 1 -p -1 -i peak/KO_merged.bed -g $FAI > peak/KO_merged.shifted.bed
bedtools shift -m 1 -p -1 -i peak/WT_merged.bed -g $FAI > peak/WT_merged.shifted.bed

bedtools genomecov -bg -strand + -5 -i peak/KO_merged.shifted.bed -g $FAI|grep '^chr' |sortBed > peak/KO_merged.shifted.plus.bedgraph &
bedtools genomecov -bg -strand - -5 -i peak/KO_merged.shifted.bed -g $FAI|grep '^chr' |sortBed > peak/KO_merged.shifted.minus.bedgraph &
bedtools genomecov -bg -strand + -5 -i peak/WT_merged.shifted.bed -g $FAI|grep '^chr' |sortBed > peak/WT_merged.shifted.plus.bedgraph &
bedtools genomecov -bg -strand - -5 -i peak/WT_merged.shifted.bed -g $FAI|grep '^chr' |sortBed > peak/WT_merged.shifted.minus.bedgraph &

./bedGraphToBigWig peak/KO_merged.shifted.plus.bedgraph $FAI peak/KO_merged.shifted.plus.bw &
./bedGraphToBigWig peak/KO_merged.shifted.minus.bedgraph $FAI peak/KO_merged.shifted.minus.bw &
./bedGraphToBigWig peak/WT_merged.shifted.plus.bedgraph $FAI peak/WT_merged.shifted.plus.bw &
./bedGraphToBigWig peak/WT_merged.shifted.minus.bedgraph $FAI peak/WT_merged.shifted.minus.bw &

~/R/3.6.1/bin/Rscript postprocess_PureCLIP_peak.R WT
~/R/3.6.1/bin/Rscript postprocess_PureCLIP_peak.R KO

~/R/3.6.1/bin/Rscript binding_sites_clean.R WT
~/R/3.6.1/bin/Rscript binding_sites_clean.R KO


# reproduce
mkdir reproduce
echo -e "file\trep\twindow\tc0\tc1\tc2\tstrand" > reproduce/crosslink_sites_reproduce_all.tsv
for F in $(find bed/*rep*sorted.bedgraph|sed 's/_rep.*//'|sort|uniq); do
    cat ${F}_rep1*minus.sorted.bedgraph|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}} END {print f"\t1\t0\t"c0"\t"c1"\t"c2"\t-"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep1*plus.sorted.bedgraph |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}} END {print f"\t1\t0\t"c0"\t"c1"\t"c2"\t+"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[23]*minus.sorted.bedgraph|bedtools window -a ${F}_rep1*minus.sorted.bedgraph -b - -w 5 -u |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t1\t5\t"c0"\t"c1"\t"c2"\t-"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[23]*plus.sorted.bedgraph |bedtools window -a ${F}_rep1*plus.sorted.bedgraph  -b - -w 5 -u |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t1\t5\t"c0"\t"c1"\t"c2"\t+"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[23]*minus.sorted.bedgraph|bedtools window -a ${F}_rep1*minus.sorted.bedgraph -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t1\t30\t"c0"\t"c1"\t"c2"\t-"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[23]*plus.sorted.bedgraph |bedtools window -a ${F}_rep1*plus.sorted.bedgraph  -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t1\t30\t"c0"\t"c1"\t"c2"\t+"}' >> reproduce/crosslink_sites_reproduce_all.tsv

    cat ${F}_rep2*minus.sorted.bedgraph|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}} END {print f"\t2\t0\t"c0"\t"c1"\t"c2"\t-"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep2*plus.sorted.bedgraph |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}} END {print f"\t2\t0\t"c0"\t"c1"\t"c2"\t+"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[13]*minus.sorted.bedgraph|bedtools window -a ${F}_rep2*minus.sorted.bedgraph -b - -w 5 -u |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t2\t5\t"c0"\t"c1"\t"c2"\t-"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[13]*plus.sorted.bedgraph |bedtools window -a ${F}_rep2*plus.sorted.bedgraph  -b - -w 5 -u |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t2\t5\t"c0"\t"c1"\t"c2"\t+"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[13]*minus.sorted.bedgraph|bedtools window -a ${F}_rep2*minus.sorted.bedgraph -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t2\t30\t"c0"\t"c1"\t"c2"\t-"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[13]*plus.sorted.bedgraph |bedtools window -a ${F}_rep2*plus.sorted.bedgraph  -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t2\t30\t"c0"\t"c1"\t"c2"\t+"}' >> reproduce/crosslink_sites_reproduce_all.tsv

    cat ${F}_rep3*minus.sorted.bedgraph|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}} END {print f"\t3\t0\t"c0"\t"c1"\t"c2"\t-"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep3*plus.sorted.bedgraph |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}} END {print f"\t3\t0\t"c0"\t"c1"\t"c2"\t+"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[12]*minus.sorted.bedgraph|bedtools window -a ${F}_rep3*minus.sorted.bedgraph -b - -w 5 -u |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t3\t5\t"c0"\t"c1"\t"c2"\t-"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[12]*plus.sorted.bedgraph |bedtools window -a ${F}_rep3*plus.sorted.bedgraph  -b - -w 5 -u |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t3\t5\t"c0"\t"c1"\t"c2"\t+"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[12]*minus.sorted.bedgraph|bedtools window -a ${F}_rep3*minus.sorted.bedgraph -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t3\t30\t"c0"\t"c1"\t"c2"\t-"}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[12]*plus.sorted.bedgraph |bedtools window -a ${F}_rep3*plus.sorted.bedgraph  -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($4>0){c0+=1}if($4>1){c1+=1}if($4>2){c2+=1}}END{print f"\t3\t30\t"c0"\t"c1"\t"c2"\t+"}' >> reproduce/crosslink_sites_reproduce_all.tsv
done

Rscript reproduced_sites.R reproduce/crosslink_sites_reproduce_all.tsv reproduce/crosslink_sites_reproduce_all.pdf

~/R/3.2.4/bin/Rscript peak_anno.R results/KO_binding_site_clean.bed /home/xfu/Gmatic7/gene/human/txdb/GRCh38_v29_txdb.sqlite figure/KO_binding_site_clean.bed_pie.pdf
~/R/3.2.4/bin/Rscript peak_anno.R results/WT_binding_site_clean.bed /home/xfu/Gmatic7/gene/human/txdb/GRCh38_v29_txdb.sqlite figure/WT_binding_site_clean.bed_pie.pdf
