#!/bin/bash

#Purpose: to process fast5 (or bam) files from nanopore run through to the end of DNAscent.
#Before you start create a folder where you would like to save run analysis, this will be </folder/to/save/run/analysis>

#Usage: bash runDNAscent.bash -f </path/to/fast5/files> -o </folder/to/save/run/analysis> [optional: -g -k -d <detect threshold> -n <output name> ] [optional: -L </path/to/bed/for/regions> | -s <INT.FRAC> ]

#optional:
# -g to do basecalling and mapping, default is off, if off requires indexed bam file called alignments.sorted, and sequencing_summary.txt to be present in -o </folder/to/save/run/analysis>. Make sure -o </folder/to/save/run/analysis> doesn't contain any files/folders that could be overwritten eg logfiles/ and bedgraphs/.
# -k to use forkSense, default off, NB forkSense has a bug, will output origin and termination bed files to pwd rather than $SAVEPATH
# -d default is 1000, the minimum for dnascent detect
# -n default is output, suggested to use if using -L or -s
# -L to generate bam for defined genomic regions and use this bam for dnascent (-L flag in samtools view)
# -s to generate subsampled bam and use this bam for dnascent (-s flag in samtools view)

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-10.1/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64

#set defualt options
BASECALL="FALSE"
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
		-o )	shift
			SAVEPATH="$1" ;;
		-g )	BASECALL="TRUE" ;;
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
	BAM="$SAVEPATH""$NAME".bam
	else
	BAM="$SAVEPATH"alignments.sorted
fi

#print variables to check

echo BASECALL = $BASECALL
echo FORKSENSE = $FORKSENSE
echo FAST5 = $FAST5
echo SAVEPATH = $SAVEPATH
echo DETECTTHRESHOLD = $DETECTTHRESHOLD
echo NAME = $NAME
echo REGION = $REGION
echo SUBSAMPLE = $SUBSAMPLE
echo "BAM to use for detect" = $BAM

#make folders and files

mkdir "$SAVEPATH"logfiles
touch "$SAVEPATH"logfiles/index_output.txt "$SAVEPATH"logfiles/detect_output.txt "$SAVEPATH"logfiles/bedgraphs_output.txt

#optional basecalling (guppy) and mapping (minimap) if starting from fast5 files
if [ "$BASECALL" != "FALSE" ]; then
	touch "$SAVEPATH"logfiles/guppy_output.txt "$SAVEPATH"logfiles/minimap_output.txt

	#use guppy to basecall fast5 files to generate fastq files, StdOut saved to guppy_ouput.txt
	/data/software_local/ont-guppy/bin/guppy_basecaller -i "$FAST5" -s "$SAVEPATH" -c /data/software_local/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg -r -x 'cuda:0' > "$SAVEPATH"logfiles/guppy_output.txt

	#tidy output
	mkdir "$SAVEPATH"fastq_files
	mv "$SAVEPATH"*.log "$SAVEPATH"logfiles/
	mv "$SAVEPATH"*.fastq "$SAVEPATH"fastq_files
	cat "$SAVEPATH"fastq_files/*.fastq > "$SAVEPATH"reads.fastq

	echo
	echo "$SAVEPATH" fastq files generated and tidied.
	echo

	#use minimap to map reads to reference, StdErr saved to minimap_ouput.txt
	/data/software_local/minimap2-2.10/minimap2 -ax map-ont -t 50 /data/workspace/rose/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna "$SAVEPATH"reads.fastq 2> "$SAVEPATH"logfiles/minimap_output.txt | samtools view -Sb - | samtools sort - -o "$SAVEPATH"alignments.sorted

	samtools index "$SAVEPATH"alignments.sorted

	echo
	echo "$SAVEPATH" reads mapped to reference.
	echo
fi

#if you want to make a smaller bam to perform DNAscent on either specific regions or a subsampke of full bam, provide arguments -L (bed file with list of regions to keep) or -s (INT.FRAC for samtools view -s subsample flag), don't use together, also provide -n <name>

if [ "$REGION" != "FALSE" ]; then
	samtools view -h -b -M -L "$REGION" -o "$SAVEPATH""$NAME".bam "$SAVEPATH"alignments.sorted
	samtools index "$SAVEPATH""$NAME".bam
	elif [ "$SUBSAMPLE" != "FALSE" ]; then
	samtools view -h -b -s "$SUBSAMPLE" -o "$SAVEPATH""$NAME".bam "$SAVEPATH"alignments.sorted
        samtools index "$SAVEPATH""$NAME".bam
fi

#run DNAscent 2.0 index and detect. StdErr saved to detect_output.txt. If you chose to make a smaller region bam then DNAscent uses this bam.
echo "DNAscent index"
/home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/bin/DNAscent index -f "$FAST5" -s "$SAVEPATH"sequencing_summary.txt -o "$SAVEPATH"index.dnascent 2> "$SAVEPATH"logfiles/index_output.txt
echo "DNAscent detect"
/home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/bin/DNAscent detect -b "$BAM" -r /data/workspace/rose/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -i "$SAVEPATH"index.dnascent -o "$SAVEPATH""$NAME".detect -t 50 --GPU 0 -l "$DETECTTHRESHOLD" 2> "$SAVEPATH"logfiles/detect_output.txt

echo
echo "$SAVEPATH" detect complete.
echo

#optional (-f) run DNAscent 2.0 forksense , StdErr saved to forkSense_output.txt
if [ "$FORKSENSE" != "FALSE" ]; then
	touch "$SAVEPATH"logfiles/forkSense_output.txt
	echo "DNAscent forkSense"
	/home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/bin/DNAscent forkSense -d "$SAVEPATH""$NAME".detect -o "$SAVEPATH""$NAME".forkSense --markOrigins --markTerminations 2> "$SAVEPATH"logfiles/forkSense_output.txt
	echo "$SAVEPATH" forksense complete.

	#Convert detect and forkSense files to bedgraphs, StdErr saved to bedgraph_output.txt
	echo "make bedgraphs"
	python /home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/utils/dnascent2bedgraph.py -d "$SAVEPATH""$NAME".detect -f "$SAVEPATH""$NAME".forkSense -o "$SAVEPATH"bedgraphs/ 2> "$SAVEPATH"logfiles/bedgraph_output.txt
	else
	#Just convert detect file to bedgraphs
	echo "make bedgraphs"
	python /home/nieduszynski/michael/development/DNAscent_v2/DNAscent_dev/utils/dnascent2bedgraph.py -d "$SAVEPATH""$NAME".detect -o "$SAVEPATH"bedgraphs/ 2> "$SAVEPATH"logfiles/bedgraph_output.txt
fi

echo
echo "$SAVEPATH" bedgraphs saved.
echo
