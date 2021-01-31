#!/bin/bash

#----------------------------------------------------------
# Copyright 2020-2021 Rosemary H C Wilson (University of Oxford) and Conrad Nieduszynski (Earlham Institute)
# Written by Rosemary H C Wilson (University of Oxford) and Conrad Nieduszynski (Earlham Institute)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

# A bash script for processing nanopore DNA sequence data through the DNAscent pipeline
# Bash script template inspiration from: https://betterdev.blog/minimal-safe-bash-script-template/
# Further bash template inspiration from: https://github.com/ralish/bash-script-template/blob/main/template.sh
# And some more advise from Google: https://google.github.io/styleguide/shellguide.html#s7-naming-conventions

# Checked with https://www.shellcheck.net

# DESC: Usage help
# ARGS: None
# OUTS: None
function usage() {
	cat << EOF

Usage: bash runDNAscent.bash -a </path/to/save/whole/run/files> -o <name_for_ouput directory> -f </path/to/fast5/files>

[ optional: -g | -m , require: -r </path/to/reference/genome>]
[ optional: -q </path/to/fastq> -k -d <detect threshold> -n <output name> -v -E -h]
[ optional: -L </path/to/bed/for/regions> | -s <INT.FRAC> ]

Purpose: to process fast5, fastq or bam files from nanopore for BrdU incorporation using DNAscent 2.0.

Before you start create a folder where you would like to save whole run files (-a below),
this will be </folder/to/save/run/files>

If using forksense run in -a directory as there was a bug (now fixed?) that origin and
termination bed files are saved to pwd then move to -o.

Can use absolute or relative paths.

Required parameters/flags:

	-a 				fastq files, sequencing summary and indexed bam of whole run saved here
	-o 				create and populate folder with any filtered indexed bam files, DNAscent
					detect and forkSense files so that you can reanalyse reads with different
					parameters without overwriting eg whole run or just specific chromosomes
	-f				fast5 files

Optional parameters/flags:

	-E|--EI				run on EI HPC. Default is to run locally (on Nieduszynski server in Oxford)
	-h|--help			Displays this help
	-v|--verbose			Displays verbose output
	-g				perform basecalling and mapping, default is off, if off requires indexed bam
					file called alignments.sorted, and sequencing_summary.txt to be present in
					-a </path/to/save/run/files> or use -m to map from fastq.
					Make sure -a doesn't contain any files/folders that could be overwritten.
	-m				to do just mapping. Default it off. If using this option it requires
					reads.fastq file in -a directory. Or chose other file with -q.
	-q				Use with -m, path to fastq file if not called reads.fastq and in -a directory.
	-r				reference genome for mapping
	-k				to use forkSense, default off
	-d				default is 1000 nts, same as default for dnascent detect
	-n				default is output, suggested to use other name especially if using -L or -s
	-L				to generate bam for defined genomic regions and use this bam for dnascent
					(-L flag in samtools view)
	-s				to generate subsampled bam and use this bam for dnascent (-s flag in samtools view)

EOF
}

# DESC: Script initialisation for use on Nieduszynski server at UoO
# ARGS: None
# OUTS: Exports locations for cuda libraries
# NOTE: This is where to add anything specific to running on Nieduszynski server.
function UoO_init() {
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-10.1/lib64
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64
	export PATH=/data/software_local/guppy_legacy/v3.6/bin:$PATH		# path to guppy
	export PATH=/data/software_local/minimap2-2.10:$PATH		# path to minimap2
	export PATH=/home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/bin:$PATH				# path to DNAscent v2
	readonly python_utils_dir="/home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/utils"	# path to DNAscent v2 utilities
	readonly guppy_model_dir="/data/software_local/guppy_legacy/v3.6/data/"										# path to guppy model files
}

# DESC: Script initialisation for use on EI HPC at NRP
# ARGS: None
# OUTS: Nothing yet
# NOTE: This is where to add anything specific to running on the EI HPC.
function EI_HPC_init() {
	source package 0e96b5e6-3f41-4d6f-91cc-1b6d7ad05ef5			# guppy - 4.0.14
	source package /tgac/software/testing/bin/minimap2-2.17		# minimap2 - 2.17
	source package 758be80b-33cc-495a-9adc-11882ab145b1			# samtools - 1.10
	source /ei/software/staging/CISSUPPORT-12154/stagingloader	# DNAscent - 2.0.2
	readonly python_utils_dir="/ei/projects/a/ac9cb897-b4c0-44d0-a54b-2ddf13310bc4/data/scripts"	# path to DNAscent v2 utilities
	readonly guppy_model_dir=""									# path to guppy model files - empty, since not required
	# Create a (hopefully) unique prefix for the names of all jobs in this 
	# particular run of the pipeline. This makes sure that runs can be
	# identified unambiguously
	run=$(uuidgen | tr '-' ' ' | awk '{print $1}')
}

# DESC: Exit script with the given message
# ARGS: $1 (required): Message to print on exit
#       $2 (optional): Exit code (defaults to 1)
# OUTS: None
function die() {
  local msg=$1
  local code=${2-1}				# default exit status 1
  echo "$msg"
  echo
  exit "$code"
}

# DESC: Parameter parser
# ARGS: $@ (optional): Arguments provided to the script
# OUTS: Variables indicating command-line parameters and options
function parse_params() {
	# default values of variables set from params
	RUNSCRIPT="UoO"
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
         		-E | --EI)
         		   	RUNSCRIPT="EI"
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
				REFGENOME="${2-}"
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
				REGION="${2-}"
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

	# check required params and not overwriting
	[[ -z "${FAST5-}" ]] && die "Missing required parameter: -f </path/to/fast5>"
	[[ -z "${RUNPATH-}" ]] && die "Missing required parameter: -a </path/to/save/whole/run/files>"
	[[ -z "${SAVEDIR-}" ]] && die "Missing required parameter: -o <name_for_ouput directory>"

	if [[ "$BASECALL" == true || "$MAPPING" == true ]] ; then
		[[ -z "${REFGENOME-}" ]] && die "Missing required parameter: -r </path/to/refgenome>"
	fi

	if [[ "$BASECALL" == true ]] ; then
		[[ -f "$RUNPATH"reads.fastq ]] && die "Reads.fastq already exists in -a directory"
	fi

	# TODO: add check that incompatible params aren't selected (e.g. -L and -s , if this really shouldn't ever be done)

	if [[ "$MAPPING" == true ]]; then
		[[ ! -f "$RUNPATH"reads.fastq && "$FASTQTEMP" == "DEFAULT" ]] \
		&& die "Require reads.fastq in -a directory or specify -q </path/to/fastq>"
		[[ -f "$RUNPATH"alignments.sorted ]] && die "Alignments.sorted already exists in -a directory"
	fi

	if [[ "$BASECALL" == "FALSE"  &&  "$MAPPING" == "FALSE" ]]; then
		[[ ! -f "$RUNPATH"alignments.sorted || ! -f "$RUNPATH"alignments.sorted.bai ]] \
		&& die "If not basecalling and/or mapping, require indexed alignments.sorted bam in -a directory"
	fi

	[[ -d "$RUNPATH""$SAVEDIR" ]] && die "Chose different name for -o save dir, already exists"

# @Conrad - how to return 0 exit code from die? This returns 1 exit code, only advice not required
	if [[ "$REGION" != "FALSE" || "$SUBSAMPLE" != "FALSE" ]]; then
		[[ "${NAME}" == "output" ]] && die "Recommended to provide -n name with -s or -L" 0
	fi

	# commented out - there aren't any required arguments
	#[[ ${#args[@]} -eq 0 ]] && die "Missing script arguments"

	return 0
}

# DESC: basecall fast5 data with guppy to generate fastq
# ARGS: None
# OUTS: None
# NOTE: This still needs to be generalised and adapted for SLURM use
function basecall_fn() {
	# check not overwriting
#	[[ -f "$RUNPATH"reads.fastq ]] && die "Exit as reads.fastq already exists in -a directory"

	local command1=(mkdir "$RUNPATH""$SAVEDIR"/logfiles/guppy_logfiles)
	local command2=(touch "$RUNPATH""$SAVEDIR"/logfiles/guppy_output.txt \
						"$RUNPATH""$SAVEDIR"/logfiles/minimap_output.txt)

	#use guppy to basecall fast5 files to generate fastq files, StdOut saved to guppy_ouput.txt
	local command3=(guppy_basecaller -i "$FAST5" -s "$RUNPATH" \
					-c "$guppy_model_dir""$guppy_model" -r \
					-x 'cuda:0' > "$RUNPATH""$SAVEDIR"/logfiles/guppy_output.txt)

	#tidy output
	local command4=(mkdir "$RUNPATH"fastq_files)
	local command5=(mv "$RUNPATH"*.log "$RUNPATH""$SAVEDIR"/logfiles/guppy_logfiles)
	local command6=(mv "$RUNPATH"*.fastq "$RUNPATH"fastq_files)
	local command7=(cat "$RUNPATH"fastq_files/*.fastq > "$RUNPATH"reads.fastq)

	if [[ "$RUNSCRIPT" == "EI" ]]; then

		# TODO: SLURM job submissions to go here

		else
		"${command1[@]}"
		"${command2[@]}"
		"${command3[@]}"
		"${command4[@]}"
		"${command5[@]}"
		"${command6[@]}"
		"${command7[@]}"		
		# check worked
		if [[ -f "$RUNPATH"reads.fastq ]]; then
			info "$RUNPATH fastq files generated and tidied"
			info
			else
			die "Exit, $RUNPATH reads.fastq not made, check guppy log files"
		fi
	fi
}

# DESC: mapping fastq data with minimap2 to generate sorted bam file (with samtools)
# ARGS: None
# OUTS: None
# NOTE: This still needs to be generalised and adapted for SLURM use
function mapping_fn() {
	#use minimap to map reads to reference, StdErr saved to minimap_ouput.txt
	local command1=(minimap2 -ax map-ont -t 50 "$REFGENOME" "$FASTQ" \
		2> "$RUNPATH""$SAVEDIR"/logfiles/minimap_output.txt \
		| samtools view -Sb - \
		| samtools sort - -o "$RUNPATH"alignments.sorted)

	local command2=(samtools index "$RUNPATH"alignments.sorted)

	if [[ "$RUNSCRIPT" == "EI" ]]; then

		# TODO: SLURM job submissions to go here

		else
		"${command1[@]}"
		"${command2[@]}"
		if [[ -f "$RUNPATH"alignments.sorted ]] ; then
			info "$FASTQ" reads mapped to reference.
			info
			else
			die "Exit, $RUNPATH alignments.sorted not made, check minimap log files"
		fi
	fi
}

# DESC: select specific regions from the bam bile with provided bed file (-L flag)
# ARGS: None
# OUTS: None
# NOTE: This still needs to be generalised and adapted for SLURM use
function regional_bam_fn() {
	samtools view -h -b -M -L "$REGION" -o "$RUNPATH""$SAVEDIR"/"$NAME".bam "$RUNPATH"alignments.sorted
	samtools index "$RUNPATH""$SAVEDIR"/"$NAME".bam
}

# DESC: subsample bam file to fraction of original size (-s flag)
# ARGS: None
# OUTS: None
# NOTE: This still needs to be generalised and adapted for SLURM use
function sub_bam_fn() {
	samtools view -h -b -s "$SUBSAMPLE" -o "$RUNPATH""$SAVEDIR"/"$NAME".bam "$RUNPATH"alignments.sorted
    samtools index "$RUNPATH""$SAVEDIR"/"$NAME".bam
}

# DESC: Function to run DNAscent 2.0 detect (and index if necessary)
# ARGS: None
# OUTS: StdErr saved to detect_output.txt.
# NOTE: If you chose to make a smaller region bam then DNAscent uses this bam.
# NOTE: This still needs to be generalised and adapted for SLURM use
function dnascent_fn() {
	if [[ ! -f "$RUNPATH"index.dnascent ]]; then
		info "DNAscent index"
		info
		DNAscent index -f "$FAST5" -s "$RUNPATH"sequencing_summary.txt -o "$RUNPATH"index.dnascent 2> "$RUNPATH""$SAVEDIR"/logfiles/index_output.txt
	fi

	info
	info "DNAscent detect"
	info
	DNAscent detect -b "$BAM" -r "$REFGENOME" -i "$RUNPATH"index.dnascent -o "$RUNPATH""$SAVEDIR"/"$NAME".detect -t 50 --GPU 0 -l "$DETECTTHRESHOLD" 2> "$RUNPATH""$SAVEDIR"/logfiles/detect_output.txt

	if [[ -f "$RUNPATH""$SAVEDIR"/"$NAME".detect ]] ; then
		info "$RUNPATH""$SAVEDIR" detect complete.
		info
		else
		die "Exit, detect file not made, check detect log files"
	fi
	info "make detect bedgraphs"
	python "$python_utils_dir"/dnascent2bedgraph.py -d "$RUNPATH""$SAVEDIR"/"$NAME".detect -o "$RUNPATH""$SAVEDIR"/"$NAME".detect.bedgraphs 2> "$RUNPATH""$SAVEDIR"/logfiles/detect_bedgraph_output.txt
}

# DESC: Function to run DNAscent 2.0 forkSense
# ARGS: None
# OUTS: StdErr saved to forkSense_output.txt
# NOTE: This still needs to be generalised and adapted for SLURM use
function forksense_fn() {
	touch "$RUNPATH""$SAVEDIR"/logfiles/forkSense_output.txt
	info "DNAscent forkSense"
	info
	DNAscent forkSense -d "$RUNPATH""$SAVEDIR"/"$NAME".detect -o "$RUNPATH""$SAVEDIR"/"$NAME".forkSense --markOrigins --markTerminations 2> "$RUNPATH""$SAVEDIR"/logfiles/forkSense_output.txt
	# to fix forksense save location bug
	mv "$RUNPATH"origins_DNAscent_forkSense.bed "$RUNPATH""$SAVEDIR"
	mv "$RUNPATH"terminations_DNAscent_forkSense.bed "$RUNPATH""$SAVEDIR"

	if [[ -f "$RUNPATH""$SAVEDIR"/"$NAME".forkSense ]]; then
		info "$RUNPATH""$SAVEDIR" forksense complete.
		info
	else
	die "Exit, forksense file not found, check forksense log file"
	fi
	info "make forksense bedgraphs"
	python "$python_utils_dir"/dnascent2bedgraph.py -f "$RUNPATH""$SAVEDIR"/"$NAME".forkSense -o "$RUNPATH""$SAVEDIR"/"$NAME".forksense.bedgraphs 2> "$RUNPATH"/"$SAVEDIR"/logfiles/forkSense_bedgraph_output.txt
}

# DESC: Generic script initialisation
# ARGS: $@ (optional): Arguments provided to the script
# OUTS: $orig_cwd: The current working directory when the script was run
#       $script_path: The full path to the script
#       $script_dir: The directory path of the script
#       $script_name: The file name of the script
#       $script_params: The original parameters provided to the script
#       $ta_none: The ANSI control code to reset all text attributes
# NOTE: $script_path only contains the path that was used to call the script
#       and will not resolve any symlinks which may be present in the path.
#       You can use a tool like realpath to obtain the "true" path. The same
#       caveat applies to both the $script_dir and $script_name variables.
function script_init() {
    # Useful paths
    readonly orig_cwd="$PWD"
    readonly script_path="${BASH_SOURCE[0]}"
    readonly script_dir="$(dirname "$script_path")"
    readonly script_name="$(basename "$script_path")"
    readonly script_params="$*"

    # Important to always set as we use it in the exit handler
    readonly ta_none="$(tput sgr0 2> /dev/null || true)"

	# Set the bam file to be used with DNAscent dependent upon whether subregion requested
	if [[ "$REGION" != "FALSE" || "$SUBSAMPLE" != "FALSE" ]]; then
		BAM="$RUNPATH""$SAVEDIR"/"$NAME".bam
		else
		BAM="$RUNPATH"alignments.sorted
	fi

	# Set the location of fastq files as specified by user
	if [[ "$FASTQTEMP" == "DEFAULT" ]]; then
		FASTQ="$RUNPATH"reads.fastq
		else
		FASTQ="$FASTQTEMP"
	fi

	#make folders and files
 	mkdir "$RUNPATH""$SAVEDIR"
	mkdir "$RUNPATH""$SAVEDIR"/logfiles
	touch "$RUNPATH""$SAVEDIR"/logfiles/index_output.txt "$RUNPATH""$SAVEDIR"/logfiles/detect_output.txt
	touch "$RUNPATH""$SAVEDIR"/logfiles/bedgraph_output.txt

    readonly guppy_model="dna_r9.4.1_450bps_fast.cfg"			# guppy model to use for basecalling
}

# DESC: Script to output information to screen (if verbose) or summary file
# ARGS: Informative text to output
# OUTS: 
# NOTE:
function info() {
	if [[ "$verbose" == true ]]; then
	    echo "INFO: $@" 
	else
		echo "INFO: $@" > "$RUNPATH""$SAVEDIR"/logfiles/summary.txt 
}

# DESC: Prints variables those variables that have been set, if verbose true
# ARGS: None
# OUTS: None
# NOTE:
function print_variables() {
	info Run script at $RUNSCRIPT
	info BASECALL and mapping = $BASECALL
	info MAPPING only = $MAPPING
	info FORKSENSE = $FORKSENSE
	info FAST5 = $FAST5
	info RUNPATH = $RUNPATH
	info SAVEDIR = $SAVEDIR
	info REFGENOME = $REFGENOME
	info DETECTTHRESHOLD = $DETECTTHRESHOLD
	info NAME = $NAME
	info REGION = $REGION
	info SUBSAMPLE = $SUBSAMPLE
	info "BAM to use for detect" = $BAM
	info "Fastq file to use" = $FASTQ
	info
}


########################################################################
# Main programme starts here
########################################################################

# First collect parameters and flags
parse_params "$@"

# Generic script initialisation that is compute independent
script_init "$@"

# Compute specific script initialisation (either UoO Nieduszynski (local) or EI HPC)
if [[ "$RUNSCRIPT" == "EI" ]]; then
	EI_HPC_init
	else
	UoO_init
fi

#print variables to check
print_variables

#optional basecalling (guppy) and mapping (minimap) if starting from fast5 files
if [[ "$BASECALL" == true ]]; then
	basecall_fn
fi

if [[ "$BASECALL" == true || "$MAPPING" == true ]]; then
	mapping_fn
fi

# if you want to make a smaller bam to perform DNAscent on either specific regions or a
# subsample of full bam, provide arguments:
# -L (bed file with list of regions to keep) or
# -s (INT.FRAC for samtools view -s subsample flag),
# don't use together, also provide -n <name>
if [[ "$REGION" != "FALSE" ]]; then
	regional_bam_fn
	elif [[ "$SUBSAMPLE" != "FALSE" ]]; then
	sub_bam_fn
fi

# run DNAscent 2.0 index (if necessary) and detect.
# StdErr saved to detect_output.txt.
# If you chose to make a smaller region bam then DNAscent uses this bam.
# Also generates bedgraphs from the detect files
dnascent_fn

if [[ -d "$RUNPATH""$SAVEDIR"/"$NAME".detect.bedgraphs ]]; then
        echo
        echo "$RUNPATH""$SAVEDIR" detect bedgraphs saved.
        else
        die "Exit, $NAME.detect.bedgraphs folder not found. Check bedgraphs logfile"
fi

# (optional) run (-f) DNAscent 2.0 forksense,
# StdErr saved to forkSense_output.txt
# Also generate bedgraphs from the forksense output
if [[ "$FORKSENSE" == true ]]; then
	forksense_fn

	if [[ -d "$RUNPATH""$SAVEDIR"/"$NAME".forksense.bedgraphs ]]; then
		info
		info "$RUNPATH""$SAVEDIR" forksense bedgraphs saved.
		else
		die "Exit, $NAME.forksense.bedgraphs folder not found. Check bedgraphs logfile"
	fi
fi

