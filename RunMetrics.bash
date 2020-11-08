#!/bin/bash

# This script generates run metrics for a run ie number of reads meeting filtering criterio and on target, also coverage in target regions.

# Useage bash RunMetrics.bash <Name of run> </path/to/bam> </path/to/region/bed/file>
# paths can be absolute or relative.

# Collect argumetns from the command line

RUN="$1"
BAM="$2"
REGIONS="$3"

# Start file for metrics and print full run read numbers
touch "$RUN".RunMetrics.txt

echo "$RUN" >> "$RUN".RunMetrics.txt
echo "Total reads for run" >> "$RUN".RunMetrics.txt
samtools view -c "$BAM" >> "$RUN".RunMetrics.txt
echo >> "$RUN".RunMetrics.txt

echo "Number of reads after secondary alignment removal and quality filtering" >> "$RUN".RunMetrics.txt
samtools view -F 0x904 -q 20 -h "$BAM" | samtools view -c >> "$RUN".RunMetrics.txt
echo >> "$RUN".RunMetrics.txt

echo "Number of on-target reads" >> "$RUN".RunMetrics.txt
samtools view -F 0x904 -q 20 -h -L "$REGIONS" "$BAM" | samtools view -c >> "$RUN".RunMetrics.txt
echo >> "$RUN".RunMetrics.txt

# for each region in bedfile of regions print read number and coverage:
while read -r CHR START END REGION; do

        echo "Region is $REGION, chrom is $CHR, start is $START, end is $END" >> "$RUN".RunMetrics.txt
	echo "Number of reads" >> "$RUN".RunMetrics.txt
	samtools view -F 0x904 -q 20 -h "$BAM" "$CHR":"$START"-"$END" | samtools view -c >> "$RUN".RunMetrics.txt
	echo >> "$RUN".RunMetrics.txt

	echo "Median coverage for $REGION" >> "$RUN".RunMetrics.txt
 	samtools depth -G 0x904 -Q 20 "$BAM" -r "$CHR":"$START"-"$END" | cut -f 3 | sort -n | awk '{ count[NR] = $1; }
	END {
  	if (NR % 2) {
	      	print count[(NR + 1) / 2];
	} else {
	       	print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0; }
	}  ' >> "$RUN".RunMetrics.txt
	echo >> "$RUN".RunMetrics.txt
done < "$REGIONS"
