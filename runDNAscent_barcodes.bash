#!/bin/bash

# Purpose: to process fast5 (or bam) files from nanopore run through to the end of DNAscent. ONLY FOR BARCODED RUNS
# Before you start create a folder where you would like to save whole run files (-a below), this will be </folder/to/save/run/files> and if you are not basecalling, a barcodes.txt file with the names of the barcode folders eg barcode01

# Usage: bash runDNAscent.bash -f </path/to/fast5/files> -a </path/to/save/whole/run/files> -o <name for ouput directory> -r </path/to/reference/genome> [ optional: -g | -m ] [ optional: -b <barcode_kit_name> -k -d <detect threshold> -n <output name> ] [optional: -L </path/to/bed/for/regions> | -s <INT.FRAC> ]

# If using forksense run in -a directory as there was a bug (now fixed?) that origin and termination bed files are saved to pwd then move to -o.
# MUST use absolute paths.

#optional:
# -g to do basecalling and mapping, default is off, if off requires indexed bam file called alignments.sorted, and sequencing_summary.txt to be present in barcode folders. Make sure -a doesn't contain any files/folders that could be overwritten.
# -m to do just mapping. Default it off. If using this option it requires "BARCODE"reads.fastq file in barcode directories.
# -b barcode option for guppy basecalling, provide barcode kit name eg "EXP-NBD104"
# -a barcode folders and sequencing summary saved here
# -o create and populate folder with any filtered indexed bam files, DNAscent detect and forkSense files so that you can reanalyse reads with different parameters without overwriting eg whole run or just specific chromosomes
# -k to use forkSense, default off
# -d default is 1000, same as default for dnascent detect
# -n default is output, suggested to use other name especially if using -L or -s
# -L to generate bam for defined genomic regions and use this bam for dnascent (-L flag in samtools view)
# -s to generate subsampled bam and use this bam for dnascent (-s flag in samtools view)

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-10.1/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64

# set paths
export PATH=/data/software_local/guppy_legacy/v3.6/bin:$PATH # path to guppy
guppy_model_dir="/data/software_local/guppy_legacy/v3.6/data/"
export PATH=/data/software_local/minimap2-2.10:$PATH            # path to minimap2
export PATH=/home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/bin:$PATH                           # path to DNAscent v2
python_utils_dir="/home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/utils"       # path to DNAscent v2 utilities

#set defualt options
BASECALL="FALSE"
BARCODE="FALSE"
MAPPING="FALSE"
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
		-b )	shift
			BARCODE="$1" ;;
		-m )	MAPPING="TRUE" ;;
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

#print variables to check

echo BASECALL and mapping = $BASECALL
echo BARCODE = $BARCODE
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

#make folders and files

mkdir "$RUNPATH""$SAVEDIR"
mkdir "$RUNPATH""$SAVEDIR"/logfiles

#optional basecalling (guppy) and mapping (minimap) if starting from fast5 files
if [ "$BASECALL" != "FALSE" ]; then
	mkdir "$RUNPATH""$SAVEDIR"/logfiles/guppy_logfiles
	touch "$RUNPATH""$SAVEDIR"/logfiles/guppy_output.txt

#use guppy to basecall fast5 files to generate fastq files, StdOut saved to guppy_ouput.txt
# barcode
	guppy_basecaller -i "$FAST5" -s "$RUNPATH" -c "$guppy_model_dir"dna_r9.4.1_450bps_fast.cfg --barcode_kits "$BARCODE" -r -x 'cuda:0' > "$RUNPATH""$SAVEDIR"/logfiles/guppy_output.txt
	#make list of barcodes to run over
	PWD_TEMP=$PWD
	cd "$RUNPATH"
	ls -d barcode* > "$RUNPATH"barcodes.txt
	cd $PWD_TEMP
	mv "$RUNPATH"*.log "$RUNPATH""$SAVEDIR"/logfiles/guppy_logfiles

	# tidy up
	for BAR in $( cat "$RUNPATH"barcodes.txt ); do
		mkdir "$RUNPATH""$BAR"/fastq_files
               	mv "$RUNPATH""$BAR"/*.fastq "$RUNPATH""$BAR"/fastq_files
               	cat "$RUNPATH""$BAR"/fastq_files/*.fastq > "$RUNPATH""$BAR"/"$BAR".reads.fastq

	done
	echo "$RUNPATH" barcode fastq files generated and tidied.
	echo

	#run DNAscent 2.0 index. StdErr saved to index_output.txt. If you chose to make a smaller region bam then DNAscent uses this bam.
	echo "DNAscent index"
	DNAscent index -f "$FAST5" -s "$RUNPATH"sequencing_summary.txt -o "$RUNPATH"index.dnascent 2> "$RUNPATH""$SAVEDIR"/logfiles/index_output.txt

fi

# for barcodes, map if required and run dnascent detect

for BAR in $( cat "$RUNPATH"barcodes.txt ); do

	# map with minimap2
	if [ "$BASECALL" != "FALSE" ] || [ "$MAPPING" != "FALSE" ]; then
		touch "$RUNPATH""$SAVEDIR"/logfiles/"$BAR".minimap_output.txt
		FASTQ="$RUNPATH""$BAR"/"$BAR".reads.fastq

		#use minimap to map reads to reference, StdErr saved to minimap_ouput.txt
		minimap2 -ax map-ont -t 50 "$REFGENOME" "$FASTQ" 2> "$RUNPATH""$SAVEDIR"/logfiles/"$BAR".minimap_output.txt | samtools view -Sb - | samtools sort - -o "$RUNPATH""$BAR"/"$BAR".alignments.sorted

		samtools index "$RUNPATH""$BAR"/"$BAR".alignments.sorted

		echo
		echo "$FASTQ" reads mapped to reference.
		echo
	fi


	#if you want to make a smaller bam to perform DNAscent on either specific regions or a subsample of full bam, provide arguments -L (bed file with list of regions to keep) or -s (INT.FRAC for samtools view -s subsample flag), don't use together, also provide -n <name>

	if [ "$REGION" != "FALSE" ]; then
		samtools view -h -b -M -L "$REGION" -o "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".bam "$RUNPATH""$BAR"/"$BAR".alignments.sorted
		samtools index "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".bam
		echo "$BAR"."$NAME".bam saved to "$RUNPATH""$SAVEDIR"
		echo
		elif [ "$SUBSAMPLE" != "FALSE" ]; then
		samtools view -h -b -s "$SUBSAMPLE" -o "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".bam "$RUNPATH""$BAR"/"$BAR".alignments.sorted
        	samtools index "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".bam
		echo "$BAR"."$NAME".bam saved to "$RUNPATH""$SAVEDIR"
		echo
	fi

	#run DNAscent 2.0 detect. StdErr saved to detect_output.txt. If you chose to make a smaller region bam then DNAscent uses this bam.

	echo "DNAscent detect"

	if [ "$REGION" != "FALSE" ] || [ "$SUBSAMPLE" != "FALSE" ]; then
        	BAM="$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".bam
        else
        	BAM="$RUNPATH""$BAR"/"$BAR".alignments.sorted
	fi

	DNAscent detect -b "$BAM" -r "$REFGENOME" -i "$RUNPATH"index.dnascent -o "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".detect -t 50 --GPU 0 -l "$DETECTTHRESHOLD" 2> "$RUNPATH""$SAVEDIR"/logfiles/"$BAR".detect_output.txt

	echo
	echo "$RUNPATH""$SAVEDIR" "$BAR" detect complete.
	echo

	#optional (-f) run DNAscent 2.0 forksense , StdErr saved to forkSense_output.txt
	if [ "$FORKSENSE" != "FALSE" ]; then
		touch "$RUNPATH""$SAVEDIR"/logfiles/"$BAR".forkSense_output.txt
		echo "DNAscent forkSense"
		echo
		DNAscent forkSense -d "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".detect -o "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".forkSense --markOrigins --markTerminations 2> "$RUNPATH""$SAVEDIR"/logfiles/"$BAR".forkSense_output.txt
		#move and rename forksense bed outputs to fix bug, may need to remove
		mv ./origins_DNAscent_forkSense.bed "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".origins_DNAscent_forkSense.bed
		mv ./terminations_DNAscent_forkSense.bed "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".terminations_DNAscent_forkSense.bed
		echo "$RUNPATH""$SAVEDIR" "$BAR" forksense complete.

		# need to rename origins and termination bed files

		#Convert detect and forkSense files to bedgraphs, StdErr saved to bedgraph_output.txt
		echo
		echo "make bedgraphs"
		python "$python_utils_dir"/dnascent2bedgraph.py -d "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".detect -f "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".forkSense -o "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".bedgraphs/ 2> "$RUNPATH"/"$SAVEDIR"/logfiles/"$BAR".bedgraph_output.txt
		else
		#Just convert detect file to bedgraphs
		echo
		echo "make bedgraphs"
		python "$python_utils_dir"/dnascent2bedgraph.py -d "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".detect -o "$RUNPATH""$SAVEDIR"/"$BAR"."$NAME".bedgraphs/ 2> "$RUNPATH""$SAVEDIR"/logfiles/"$BAR".bedgraph_output.txt
	fi

	echo
	echo "$RUNPATH""$SAVEDIR" "$BAR" bedgraphs saved.
	echo
done
