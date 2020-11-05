#!/bin/bash
#SBATCH --job-name=BWA_MEMmap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1
#SBATCH --time=7-00:00:00
#SBATCH --output=/data/nieduszynski/NGS/2020_07_08_JC_RW_ILL_HeLa_BrdU_lowchases/bwa/bwa_out.%j
#SBATCH --error=/data/nieduszynski/NGS/2020_07_08_JC_RW_ILL_HeLa_BrdU_lowchases/bwa/bwa_err.%j

#usage sbatch bwa_map.bash <samples.txt> NB <samples.txt> should have stem of the filename with no extension

#should be used with samtools version 1.10, module load SAMTOOLS/1.10
#if default samtools already loaded, module unload SAMTOOLS, module load SAMTOOLS/1.10
#to use with older versions of samtools remove -M flag in 2nd samtools view command, and for samtools markdups command replace -f flag (and associated filename to save stats to) with -s, stats will be printed to screen.

#For use on barcode trimmed file in trim/ directory wihin run directory, run script in bwa/ directory under run directory

for file in $(cat $1)
do
echo $file

#map trimmed reads to specified human reference genome (no alt version) using bwa-mem, sort and save as a bam file
bwa mem -M -t 16 /data/nieduszynski/RESOURCES/REFERENCE_GENOMES/STAR/ncbi-noALT-2020_03_17/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ../trim/$file.fastq.gz | samtools sort -@ 15 -o $file.bwa.sorted.bam -
echo $file.bwa.sorted.bam saved
echo ""

#index the bam file, print mapping stats to <sample>.mappingStats.txt file
samtools index $file.bwa.sorted.bam
touch $file.mappingStats.txt
samtools flagstat $file.bwa.sorted.bam > $file.mappingStats.txt
echo "uniquely mapped reads " >> $file.mappingStats.txt
samtools view -h -@ 15 -F 3844 -q 1 $file.bwa.sorted.bam | grep -v -E "SA:Z:|XA:Z:" | wc -l >> $file.mappingStats.txt

#select uniquely mapping reads, mark and remove duplicaes (duplicate stats saved in <sample>.markdups.txt file, generate genome coverage for the 5' end of reads, keep only bases with >= 1 read mapping there and save as a bed file.
samtools view -h -@ 15 -F 3844 -q 1 -M -L /data/nieduszynski/NGS/2020_03_04_ILL_RW_HeLa_Chases/genomeWindows/hg38_10kbWindows.bed $file.bwa.sorted.bam | grep -v -E "SA:Z:|XA:Z:" | samtools view -@ 15 -b - | samtools markdup -r -f $file.markdups.txt - - | bedtools genomecov -5 -d -ibam stdin | awk 'BEGIN {OFS="\t"} {if ($3>0) print $1,$2,$2,"name",$3}' > $file.coverage.bed
echo $file.coverage.bed saved
done
