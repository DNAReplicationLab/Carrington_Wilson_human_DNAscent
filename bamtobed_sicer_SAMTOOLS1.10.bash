#!/bin/bash

# usage bash bamtobed_sicer.bash  <samples.txt> <path/to/files> <path/to/chr/bed/file> <path/to/blacklist/bed/file>

# NB <samples.txt> should have stem of the filename with no extension

PATH=$2
CHRBED=$3
BLACKLIST=$4

echo SAMPLES = $1
echo PATH = $PATH
echo CHRBED = $CHRBED
echo BLACKLIST = $BLACKLIST
echo

#should be used with samtools version 1.10, module load SAMTOOLS/1.10
#if default samtools already loaded, module unload SAMTOOLS, module load SAMTOOLS/1.10
#to use with older versions of samtools remove -M flag in 1st samtools view command

#For use on .bwa.sorted.bam files. Run in directory to save files.

for file in $( cat "$1" );
	do
	echo $file
	#filter bam for only uniquely mapped reads and save as uniq.bam NB not dedup'd but sicer does this
	samtools view -h -@ 15 -F 3844 -q 1 -M -L "$CHRBED" "$PATH""$file".bwa.sorted.bam | grep -v -E "SA:Z:|XA:Z:" | samtools view -@ 15 -b -h - > $file.uniq.bam
	echo $file ".uniq.bam saved"
	#remove areas in blacklist, sort and save as .bed
	bedtools bamtobed -i $file.uniq.bam | bedtools intersect -a stdin -b "$BLACKLIST" -v | sort -k 1,1 -k2,2n > $file.uniq.bed
	echo $file ".uniq.bed saved"
done
