#!/bin/bash

# usage bash IPoverINPUT.bash -o </path/to/folder/with/files> -j </path/to/forjoin.txt> -w <window_name> -c </path/to/chrom/sizes.bed>

# Script takes intermediate windowed coverage bed files with extra column for joining, output from genomecoverageToBigWig.bash and generates IPoverINPUT bedfiles and BigWig files
# NB <forjoin.txt> should have 3x tab separated columns with stem (before first .) of IP filesname, stem of INPUT filename and sample name
# eg IP_0P	INPUT_OP	OP

# Collect command line arguments

while [ "$1" != "" ]; do
	case $1 in
		-o )	shift
			SAVEPATH="$1" ;;
		-j )	shift
			FORJOIN="$1" ;;
		-w )	shift
			WINDOW="$1" ;;
		-c )	shift
			CHROMSIZE="$1" ;;
	esac
		shift
done

echo SAVEPATH = $SAVEPATH
echo FORJOIN = $FORJOIN
echo WINDOW = $WINDOW
echo CHROMSIZE = $CHROMSIZE
echo

echo "Generate IPoverINPUT bed files and bigwig files for window size $WINDOW samples in $FORJOIN and save at $SAVEPATH using $CHROMSIZE."
echo

# For the IP and INPUT file pairs listed in $FORJOIN file, generate IPoverINPUT files
while read -r IP INPUT NAME; do
	echo "$IP, $INPUT, $NAME"
	# make temp bed file with IP and INPUT files joined together
        join -1 5 -2 5 "$SAVEPATH""$IP"."$WINDOW".forjoin "$SAVEPATH""$INPUT"."$WINDOW".forjoin > "$SAVEPATH""$NAME".IPoverINPUT."$WINDOW".joined
	# Make new bed file with IPoverINPUT column
	awk '{print $2, $3, $4, $5/$9}' "$SAVEPATH""$NAME".IPoverINPUT."$WINDOW".joined | sort -k 1,1 -k2,2n > "$SAVEPATH""$NAME".IPoverINPUT."$WINDOW".bed
	# save bigwig
	#bedGraphToBigWig "$SAVEPATH""$NAME".IPoverINPUT."$WINDOW".bed "$CHROMSIZE" "$SAVEPATH""$NAME".IPoverINPUT."$WINDOW".bw
	#echo "$NAME joined bedfile, IPoverINPUT bedfile and BigWig saved"
	echo
done < "$FORJOIN"

