#!/bin/bash

# usage bash map_SICERpeaks.bash <my_files.txt> <path/to/files> <path/to/blacklist> <path/to/sicer/peaks/bed>

# NB my_files.txt should have stem of the filename with no extension
# use with .coverage.bed files

#should be used with samtools version 1.10, module load SAMTOOLS/1.10

#remove genome positions that overlap with blacklist, count reads in SICER2 w2000f400 significant peak call regions

for file in $( cat $1 )
	do
	echo $file
	bedtools intersect -a $2$file.coverage.bed -b $3 -v | bedtools map -o sum -null 0 -a $4 -b stdin > $file.SICERpeaks.bed
	echo $file ".SICERpeaks.bed saved"
done
