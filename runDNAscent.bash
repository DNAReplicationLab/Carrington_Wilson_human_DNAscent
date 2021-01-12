#!/bin/bash

#----------------------------------------------------------
# Copyright 2020-2021 Rose Wilson (University of Oxford) and Conrad Nieduszynski (Earlham Institute)
# Written by Rose Wilson (University of Oxford) and Conrad Nieduszynski (Earlham Institute)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

# A bash script for processing nanopore DNA sequence data through the DNAscent pipeline
# Bash script template inspiration from: https://betterdev.blog/minimal-safe-bash-script-template/
# Further bash template inspiration from: https://github.com/ralish/bash-script-template/blob/main/template.sh

# DESC: Usage help
# ARGS: None
# OUTS: None
function usage() {
	cat << EOF
Usage: bash runDNAscent.bash -f </path/to/fast5/files> -a </path/to/save/whole/run/files>
-o <name for ouput directory> -r </path/to/reference/genome> [ optional: -g | -m ]
[ optional: -q </path/to/fastq> -k -d <detect threshold> -n <output name> ]
[optional: -L </path/to/bed/for/regions> | -s <INT.FRAC> ]

Purpose: to process fast5 (or bam) files from nanopore run through to the end of DNAscent.

Before you start create a folder where you would like to save whole run files (-a below),
this will be </folder/to/save/run/files>

If using forksense run in -a directory as there was a bug (now fixed?) that origin and
termination bed files are saved to pwd then move to -o.

Can use absolute or relative paths.

Available options:

-g			to do basecalling and mapping, default is off, if off requires indexed bam
			file called alignments.sorted, and sequencing_summary.txt to be present in
			-a </path/to/save/run/files>. Make sure -a doesn't contain any files/folders
			that could be overwritten.
-m			to do just mapping. Default it off. If using this option it requires
			reads.fastq file in -a directory. Or chose other file with -q.
-q			Use with -m, path to fastq file if not called reads.fastq and in -a directory.
-a			fastq files, sequencing summary and indexed bam of whole run saved here
-o			create and populate folder with any filtered indexed bam files, DNAscent
			detect and forkSense files so that you can reanalyse reads with different
			parameters without overwriting eg whole run or just specific chromosomes
-k			to use forkSense, default off
-d			default is 1000 nts, same as default for dnascent detect
-n			default is output, suggested to use other name especially if using -L or -s
-L			to generate bam for defined genomic regions and use this bam for dnascent
			(-L flag in samtools view)
-s			to generate subsampled bam and use this bam for dnascent (-s flag in samtools view)
EOF
}

# DESC: Exit script with the given message
# ARGS: $1 (required): Message to print on exit
#       $2 (optional): Exit code (defaults to 1)
# OUTS: None
function die() {
  local msg=$1
  local code=${2-1}				# default exit status 1
  echo "$msg"
  exit "$code"
}

# DESC: Parameter parser
# ARGS: $@ (optional): Arguments provided to the script
# OUTS: Variables indicating command-line parameters and options
function parse_params() {
	# default values of variables set from params
	FAST5=''
	BASECALL="FALSE"
	MAPPING="FALSE"
	FASTQTEMP="DEFAULT"
	FORKSENSE="FALSE"
	DETECTTHRESHOLD=1000
	NAME="output"
	REGION="FALSE"
	SUBSAMPLE="FALSE"

	#collect command line arguments
	while :; do
		case "${1-}" in
			-h | --help)
				usage
				exit 0
				;;
            -v | --verbose)
                verbose=true
                ;;
			-f )	
				FAST5="${2-}"
				shift
				;;
			-a )
				RUNPATH="${2-}"
				shift
				;;
			-o )
				SAVEDIR="${2-}"
				shift
				;;
			-g )
				BASECALL=true
				;;
			-m )
				MAPPING=true
				;;
			-q )
				FASTQTEMP="${2-}"
				shift
				;;
			-r )
				REFGENOME"${2-}"
				shift
				;;
			-k )
				FORKSENSE=true
				;;
			-d )
				DETECTTHRESHOLD="${2-}"
				shift
				;;
			-n )
				NAME="${2-}"
				shift
				;;
			-L )
				REGION="$1""${2-}"
				shift
				;;
			-s )
				SUBSAMPLE="${2-}"
				shift
				;;
			-?*)
				die "Unknown option: $1"
				;;
			*)
				break
				;;
		esac
		shift
	done
	
	args=("$@")

	# check required params and arguments
	[[ -z "${FAST5-}" ]] && die "Missing required parameter: fast5"
	[[ ${#args[@]} -eq 0 ]] && die "Missing script arguments"
	
	return 0
}

parse_params "$@"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-10.1/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64


if [ "$REGION" = true ] || [ "$SUBSAMPLE" = true ]; then
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
if [ "$BASECALL" = true ]; then
	mkdir "$RUNPATH""$SAVEDIR"/logfiles/guppy_logfiles
	touch "$RUNPATH""$SAVEDIR"/logfiles/guppy_output.txt "$RUNPATH""$SAVEDIR"/logfiles/minimap_output.txt

	#use guppy to basecall fast5 files to generate fastq files, StdOut saved to guppy_ouput.txt
	/data/software_local/ont-guppy/bin/guppy_basecaller -i "$FAST5" -s "$RUNPATH" -c /data/software_local/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg -r -x 'cuda:0' > "$RUNPATH""$SAVEDIR"/logfiles/guppy_output.txt

	#tidy output
	mkdir "$RUNPATH"fastq_files
	mv "$RUNPATH"*.log "$RUNPATH""$SAVEDIR"/logfiles/guppy_logfiles
	mv "$RUNPATH"*.fastq "$RUNPATH"fastq_files
	cat "$RUNPATH"fastq_files/*.fastq > "$RUNPATH"reads.fastq

	echo
	echo "$RUNPATH" fastq files generated and tidied.
	echo
fi

if [ "$BASECALL" = true ] || [ "$MAPPING" = true ]; then
	#use minimap to map reads to reference, StdErr saved to minimap_ouput.txt
	/data/software_local/minimap2-2.10/minimap2 -ax map-ont -t 50 "$REFGENOME" "$FASTQ" 2> "$RUNPATH""$SAVEDIR"/logfiles/minimap_output.txt | samtools view -Sb - | samtools sort - -o "$RUNPATH"alignments.sorted

	samtools index "$RUNPATH"alignments.sorted

	echo
	echo "$FASTQ" reads mapped to reference.
	echo
fi

#if you want to make a smaller bam to perform DNAscent on either specific regions or a subsample of full bam, provide arguments -L (bed file with list of regions to keep) or -s (INT.FRAC for samtools view -s subsample flag), don't use together, also provide -n <name>

if [ "$REGION" = true ]; then
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
if [ "$FORKSENSE" = true ]; then
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
