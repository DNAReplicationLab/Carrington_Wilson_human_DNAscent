#!/bin/bash

# Generates intermediate BigWig files for coverage of IP or INPUT files for visualisation and generates intermediate file for performing IPoverINPUT step, both from genomecoverage bed files generated from bwa_map.bash..

# usage bash gencoverageToBigWig.bash -o </path/to/folder/with/files/and/for/save> -s </path/to/my_samples.txt> -b </path/to/blacklist> -w </path/to/windows/bed/file> -n <window_name> -c </path/to/chromosome/size/bed/file>

# NB my_samples.txt should have stem of the filename with no extension
# do not use with sbatch - if need to,  need to change sort -k1,1 -k2,2n to samtools sort.

# collect command line varibles
while [ "$1" != "" ]; do
        case $1 in
                -o )    shift
                        SAVEPATH="$1" ;;
                -s )    shift
                        SAMPLES="$1" ;;
                -b )    shift
                        BLACKLIST="$1" ;;
                -w )    shift
                        WINDOWS="$1" ;;
                -n )    shift
                        WINDOW="$1" ;;
	        -c )	shift
			CHROMSIZE="$1" ;;
esac
                shift
done

echo SAVEPATH = $SAVEPATH
echo SAMPLES = $SAMPLES
echo BLACKLIST = $BLACKLIST
echo WINDOWS = $WINDOWS
echo WINDOWNAME = $WINDOW
echo CHROMSIZE = $CHROMSIZE

echo "Generate bigwig files for samples in $SAMPLES and save at $SAVEPATH. Remove $BLACKLIST and map into $WINDOW size windows using $WINDOWSand $CHROMSIZE."

# Porcess files listed in $SAMPLES
for SAMPLE in $(cat "$SAMPLES")
	do
	echo "$SAMPLE"
	# remove any positions in blacklist, map into windows of chosen size, sort, awk command removes any empty windows
	bedtools intersect -a "$SAVEPATH""$SAMPLE".coverage.bed -b "$BLACKLIST" -v | bedtools map -o sum -null 0 -a "$WINDOWS" -b stdin | sort -k1,1 -k2,2n | awk 'BEGIN {OFS="\t"} {if ($4>0) print $1,$2,$3,$4}' > "$SAVEPATH""$SAMPLE".sorted."$WINDOW".bed
	# convert to bigwig files for viewing in IGV
	bedGraphToBigWig "$SAVEPATH""$SAMPLE".sorted."$WINDOW".bed "$CHROMSIZE" "$SAVEPATH""$SAMPLE"."$WINDOW".bw
	# Generate intermediate windowed bed file for all samples ready for IPoverINPUT step (making join_by column)
        awk '{print $1, $2, $3, $4, $1"_"$2}' "$SAVEPATH""$SAMPLE".sorted."$WINDOW".bed | sort -k 5 >  "$SAVEPATH""$SAMPLE"."$WINDOW".forjoin
	echo "$SAMPLE.$WINDOW.bw saved"
	echo
done
