#!/bin/bash

# Purpose: to process fast5 (or bam) files from nanopore run through to the end of DNAscent.
# Before you start create a folder where you would like to save whole run files (-a below), this will be </folder/to/save/run/files>

# Usage: bash runDNAscent.bash -f </path/to/fast5/files> -a </path/to/save/whole/run/files> -o <name for ouput directory> -r </path/to/reference/genome> [ optional: -g | -m ] [ optional: -q </path/to/fastq> -k -d <detect threshold> -n <output name> ] [optional: -L </path/to/bed/for/regions> | -s <INT.FRAC> ]

# If using forksense run in -a directory as there was a bug (now fixed?) that origin and termination bed files are saved to pwd then move to -o.
# Can use absolute or relative paths.

#optional:
# -g to do basecalling and mapping, default is off, if off requires indexed bam file called alignments.sorted, and sequencing_summary.txt to be present in -a </path/to/save/run/files>. Make sure -a doesn't contain any files/folders that could be overwritten.
# -m to do just mapping. Default it off. If using this option it requires reads.fastq file in -a directory. Or chose other file with -q.
# -q Use with -m, path to fastq file if not called reads.fastq and in -a directory.
# -a fastq files, sequencing summary and indexed bam of whole run saved here
# -o create and populate folder with any filtered indexed bam files, DNAscent detect and forkSense files so that you can reanalyse reads with different parameters without overwriting eg whole run or just specific chromosomes
# -k to use forkSense, default off
# -d default is 1000, same as default for dnascent detect
# -n default is output, suggested to use other name especially if using -L or -s
# -L to generate bam for defined genomic regions and use this bam for dnascent (-L flag in samtools view)
# -s to generate subsampled bam and use this bam for dnascent (-s flag in samtools view)

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-10.1/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64

# path to guppy legacy as current guppy not working with GPU currently
export PATH=/data/software_local/guppy_legacy/v3.6/bin:$PATH
guppy_model_dir="/data/software_local/guppy_legacy/v3.6/data/"

#set defualt options
BASECALL="FALSE"
MAPPING="FALSE"
FASTQTEMP="DEFAULT"
FORKSENSE="FALSE"
DETECTTHRESHOLD=1000
NAME="output"
REGION="FALSE"
SUBSAMPLE="FALSE"

#collect command line arguments
while [ "$1" != "" ]; do
	case $1 in
		-f )	shift
			FAST5="$1" ;;
		-a )	shift
			RUNPATH="$1" ;;
		-o )	shift
			SAVEDIR="$1" ;;
		-g )	BASECALL="TRUE" ;;
		-m )	MAPPING="TRUE" ;;
		-q )	shift
			FASTQTEMP="$1" ;;
		-r )	shift
			REFGENOME="$1" ;;
		-k )	FORKSENSE="TRUE" ;;
		-d )	shift
			DETECTTHRESHOLD="$1" ;;
		-n )	shift
			NAME="$1" ;;
		-L )	shift
			REGION="$1" ;;
		-s )	shift
			SUBSAMPLE="$1" ;;
	esac
	shift
done


if [ "$REGION" != "FALSE" ] || [ "$SUBSAMPLE" != "FALSE" ]; then
	BAM="$RUNPATH""$SAVEDIR"/"$NAME".bam
	else
	BAM="$RUNPATH"alignments.sorted
fi

if [ "$FASTQTEMP" == "DEFAULT" ]; then
	FASTQ="$RUNPATH"reads.fastq
	else
	FASTQ="$FASTQTEMP"
fi

#print variables to check

echo BASECALL and mapping = $BASECALL
echo MAPPING only = $MAPPING
echo FORKSENSE = $FORKSENSE
echo FAST5 = $FAST5
echo RUNPATH = $RUNPATH
echo SAVEDIR = $SAVEDIR
echo REFGENOME = $REFGENOME
echo DETECTTHRESHOLD = $DETECTTHRESHOLD
echo NAME = $NAME
echo REGION = $REGION
echo SUBSAMPLE = $SUBSAMPLE
echo "BAM to use for detect" = $BAM
echo "Fastq file to use" = $FASTQ

#make folders and files

mkdir "$RUNPATH""$SAVEDIR"
mkdir "$RUNPATH""$SAVEDIR"/logfiles
touch "$RUNPATH""$SAVEDIR"/logfiles/index_output.txt "$RUNPATH""$SAVEDIR"/logfiles/detect_output.txt "$RUNPATH""$SAVEDIR"/logfiles/bedgraph_output.txt

#optional basecalling (guppy) and mapping (minimap) if starting from fast5 files
if [ "$BASECALL" != "FALSE" ]; then
	mkdir "$RUNPATH""$SAVEDIR"/logfiles/guppy_logfiles
	touch "$RUNPATH""$SAVEDIR"/logfiles/guppy_output.txt "$RUNPATH""$SAVEDIR"/logfiles/minimap_output.txt

	#use guppy to basecall fast5 files to generate fastq files, StdOut saved to guppy_ouput.txt
	guppy_basecaller -i "$FAST5" -s "$RUNPATH" -c "$guppy_model_dir"dna_r9.4.1_450bps_fast.cfg -r -x 'cuda:0' > "$RUNPATH""$SAVEDIR"/logfiles/guppy_output.txt

	#tidy output
	mkdir "$RUNPATH"fastq_files
	mv "$RUNPATH"*.log "$RUNPATH""$SAVEDIR"/logfiles/guppy_logfiles
	mv "$RUNPATH"*.fastq "$RUNPATH"fastq_files
	cat "$RUNPATH"fastq_files/*.fastq > "$RUNPATH"reads.fastq

	echo
	echo "$RUNPATH" fastq files generated and tidied.
	echo
fi

if [ "$BASECALL" != "FALSE" ] || [ "$MAPPING" != "FALSE" ]; then
	#use minimap to map reads to reference, StdErr saved to minimap_ouput.txt
	/data/software_local/minimap2-2.10/minimap2 -ax map-ont -t 50 "$REFGENOME" "$FASTQ" 2> "$RUNPATH""$SAVEDIR"/logfiles/minimap_output.txt | samtools view -Sb - | samtools sort - -o "$RUNPATH"alignments.sorted

	samtools index "$RUNPATH"alignments.sorted

	echo
	echo "$FASTQ" reads mapped to reference.
	echo
fi

#if you want to make a smaller bam to perform DNAscent on either specific regions or a subsample of full bam, provide arguments -L (bed file with list of regions to keep) or -s (INT.FRAC for samtools view -s subsample flag), don't use together, also provide -n <name>

if [ "$REGION" != "FALSE" ]; then
	samtools view -h -b -M -L "$REGION" -o "$RUNPATH""$SAVEDIR"/"$NAME".bam "$RUNPATH"alignments.sorted
	samtools index "$RUNPATH""$SAVEDIR"/"$NAME".bam
	elif [ "$SUBSAMPLE" != "FALSE" ]; then
	samtools view -h -b -s "$SUBSAMPLE" -o "$RUNPATH""$SAVEDIR"/"$NAME".bam "$RUNPATH"alignments.sorted
        samtools index "$RUNPATH""$SAVEDIR"/"$NAME".bam
fi

#run DNAscent 2.0 index and detect. StdErr saved to detect_output.txt. If you chose to make a smaller region bam then DNAscent uses this bam.
echo "DNAscent index"
/home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/bin/DNAscent index -f "$FAST5" -s "$RUNPATH"sequencing_summary.txt -o "$RUNPATH"index.dnascent 2> "$RUNPATH""$SAVEDIR"/logfiles/index_output.txt
echo "DNAscent detect"
/home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/bin/DNAscent detect -b "$BAM" -r "$REFGENOME" -i "$RUNPATH"index.dnascent -o "$RUNPATH""$SAVEDIR"/"$NAME".detect -t 50 --GPU 0 -l "$DETECTTHRESHOLD" 2> "$RUNPATH""$SAVEDIR"/logfiles/detect_output.txt

echo
echo "$RUNPATH""$SAVEDIR" detect complete.
echo

#optional (-f) run DNAscent 2.0 forksense , StdErr saved to forkSense_output.txt
if [ "$FORKSENSE" != "FALSE" ]; then
	touch "$RUNPATH""$SAVEDIR"/logfiles/forkSense_output.txt
	echo "DNAscent forkSense"
	/home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/bin/DNAscent forkSense -d "$RUNPATH""$SAVEDIR"/"$NAME".detect -o "$RUNPATH""$SAVEDIR"/"$NAME".forkSense --markOrigins --markTerminations 2> "$RUNPATH""$SAVEDIR"/logfiles/forkSense_output.txt
	echo "$RUNPATH""$SAVEDIR" forksense complete.

	#Convert detect and forkSense files to bedgraphs, StdErr saved to bedgraph_output.txt
	echo "make bedgraphs"
	python /home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/utils/dnascent2bedgraph.py -d "$RUNPATH""$SAVEDIR"/"$NAME".detect -f "$RUNPATH""$SAVEDIR"/"$NAME".forkSense -o "$RUNPATH""$SAVEDIR"/bedgraphs/ 2> "$RUNPATH"/"$SAVEDIR"/logfiles/bedgraph_output.txt
	else
	#Just convert detect file to bedgraphs
	echo "make bedgraphs"
	python /home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/utils/dnascent2bedgraph.py -d "$RUNPATH""$SAVEDIR"/"$NAME".detect -o "$RUNPATH""$SAVEDIR"/bedgraphs/ 2> "$RUNPATH""$SAVEDIR"/logfiles/bedgraph_output.txt
fi

echo
echo "$RUNPATH""$SAVEDIR" bedgraphs saved.
echo
