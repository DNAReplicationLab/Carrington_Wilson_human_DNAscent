#!/bin/bash

#SBATCH --job-name=BWA_map
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1
#SBATCH --time=7-00:00:00
#SBATCH --output=./bwa_map_out.%j
#SBATCH --error=./bwa_map_err.%j

# usage sbatch bwa_map.bash -f </path/to/folder/with/fastq/files> -s </path/to/samples.txt> -o </path/to/save/files/> -g </path/to/reference/genome> -L </path/to/bed/of/of/genome/region/to/include>

# NB <samples.txt> should have stem of the filename with no extension
# For use on barcode trimmed fastq.gz files
# Slurm output files will be saved in pwd

# this section currently doesn't work with SLURM:
# so need to manually load BWA, module unload SAMTOOLS, load SAMTOOLS/1.10 and BEDTOOLS

# load required modules and set variables
#module load SLURM/5.08.6
#module load BWA/0.7.12 && BWA="bwa"
# requires SAMTOOLS/1.10 to use -M flag for samtools view. To use with older versions, remove -M flag from samtools view command,
# and for samtools markdups command replace -f flag (and associated filename to save stats to) with -s, stats will be printed to screen.
#module unload SAMTOOLS
#module load SAMTOOLS/1.10 && SAMTOOLS="samtools"
#module load PICARDTOOLS/6.15 && PICARDTOOLS="/software/BIOINFORMATICS/PICARDTOOLS/picard-06.15/dist/picard.jar"
#module load BEDTOOLS/2.26.0 && BEDTOOLS="bedtools"
#echo
#echo "Required modules are loaded"
#echo

# collect command line varibles
while [ "$1" != "" ]; do
	case $1 in
		-f )    shift
                        FASTQ="$1" ;;
                -s )    shift
                        SAMPLES="$1" ;;
		-o )	shift
			SAVEPATH="$1" ;;
		-g )	shift
			GENOME="$1" ;;
		-L )	shift
			INCLUDE="$1" ;;
        esac
                shift
done

echo "Map reads for samples in $SAMPLES FASTQ files in $FASTQ and save at $SAVEPATH. Use $GENOME reference genome and $INCLUDE."
echo

for SAMPLE in $(cat "$SAMPLES")
do
echo "$SAMPLE"
echo

#map trimmed reads to specified human reference genome (no alt version) using bwa-mem, sort and save as a bam file
"$BWA" mem -M -t 16 "$GENOME" "$FASTQ""$SAMPLE".trim.fastq.gz | "$SAMTOOLS" sort -@ 15 -o "$SAVEPATH""$SAMPLE".bwa.sorted.bam -
echo "$SAMPLE".bwa.sorted.bam saved
echo

#index the bam file, print mapping stats to <$SAMPLE>.mappingStats.txt file
"$SAMTOOLS" index "$SAVEPATH""$SAMPLE".bwa.sorted.bam
touch "$SAVEPATH""$SAMPLE".mappingStats.txt
"$SAMTOOLS" flagstat "$SAVEPATH""$SAMPLE".bwa.sorted.bam > "$SAVEPATH""$SAMPLE".mappingStats.txt
echo "uniquely mapped reads " >> "$SAVEPATH""$SAMPLE".mappingStats.txt
"$SAMTOOLS" view -h -@ 15 -F 3844 -q 1 "$SAVEPATH""$SAMPLE".bwa.sorted.bam | grep -v -E "SA:Z:|XA:Z:" | wc -l >> "$SAVEPATH""$SAMPLE".mappingStats.txt

#select uniquely mapping reads, mark and remove duplicaes (duplicate stats saved in <sample>.markdups.txt file, generate genome coverage for the 5' end of reads, keep only bases with >= 1 read mapping there and save as a bed file.
"$SAMTOOLS" view -h -@ 15 -F 3844 -q 1 -M -L "$INCLUDE" "$SAVEPATH""$SAMPLE".bwa.sorted.bam | grep -v -E "SA:Z:|XA:Z:" | "$SAMTOOLS" view -@ 15 -b - | "$SAMTOOLS" markdup -r -f "$SAVEPATH""$SAMPLE".markdups.txt - - | "$BEDTOOLS" genomecov -5 -d -ibam stdin | awk 'BEGIN {OFS="\t"} {if ($3>0) print $1,$2,$2,"name",$3}' > "$SAVEPATH""$SAMPLE".coverage.bed
echo "$SAMPLE".coverage.bed saved
echo
done
