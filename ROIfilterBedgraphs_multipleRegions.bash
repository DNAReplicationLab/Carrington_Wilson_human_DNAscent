#!/bin/bash

# Takes detect file and bedgraphs from dnascent and copies bedgraphs corresponding to region of interest to a new folder

# usage bash ROIfilterBedgraphs.bash -d </path/to/.detect/file> -b </path/to/bedgraphs> -s </save/path> -r </path/to/bed/for/regions> optional [-w <value of ROI extension> ]]
# -r bedfile in form chrom start end regionName for your regions of interest
# optional -w allows specifying of optional value to enxent ROI at each side, give value in base pairs eg for 100kb -w 100000

# NB paths MUST be full paths

# variables
WIDE=0
WD="$PWD"

# collect command line varibles
while [ "$1" != "" ]; do
        case $1 in
                -d )    shift
                        DETECT="$1" ;;
		-b )	shift
			BEDGRAPHS="$1" ;;
		-s )	shift
			SAVEPATH="$1" ;;
		-r )	shift
			REGIONFILE="$1" ;;
		-w )	shift
			WIDE="$1" ;;
	esac
		shift
done

# print variables
echo
echo "Detect bedgraphs to be copied to $SAVEPATH for regions in $REGIONFILE from $BEDGRAPHS using $DETECT"

# extend region if -w
if [ "$WIDE" != 0 ]; then
	echo
	echo "Region will be extended by $WIDE each side."
fi

# copy all bedgraphs to single folder
mkdir "$SAVEPATH"allBedgraphs
cp "$BEDGRAPHS"*/*.bedgraph "$SAVEPATH"allBedgraphs/

# for each region in bedfile of regions:
while read -r CHR START END REGION; do

	# extend region if -w
	if [ "$WIDE" != 0 ]; then
        	START=$(expr "$START" - "$WIDE")
        	END=$(expr "$END" + "$WIDE")
	fi
        echo
	echo "Region is $REGION, chrom is $CHR, start is $START, end is $END"
	echo
        mkdir "$SAVEPATH""$REGION".bedgraphs

	# make list of reads that fall within ROI + 100kb either side
	mkdir "$SAVEPATH""$REGION"temp
	grep "^>" "$DETECT" > "$SAVEPATH""$REGION"temp/readNamesDetect.txt
	cat "$SAVEPATH""$REGION"temp/readNamesDetect.txt | grep "$CHROM" > "$SAVEPATH""$REGION"temp/readNamesChr.txt

	#pass start and end coordinates to awk, keep lines where reads are within ROI genome coordinates.
	awk -v start="$START" -v end="$END" 'start < $4 && $3 < end' "$SAVEPATH""$REGION"temp/readNamesChr.txt > "$SAVEPATH""$REGION".readNames.txt
	cat "$SAVEPATH""$REGION".readNames.txt | cut -f2 -d ">" | cut -f1 -d " " > "$SAVEPATH""$REGION"temp/readNamesROI.txt
	echo
	echo "Reads identified for $REGION, $CHR, start $START, end $END"

	#copy ROI bedgraphs to ROI folder
	cd "$SAVEPATH"allBedgraphs/
	cp $(ls . | grep -f "$SAVEPATH""$REGION"temp/readNamesROI.txt) "$SAVEPATH""$REGION".bedgraphs/
	cd "$WD"

	#remove temp folder, comment out line below if want to check intermediate files
	rm -r "$SAVEPATH""$REGION"temp/

	echo
	echo "Reads copied for $REGION, $CHR, start $START, end $END"
done < "$REGIONFILE"

# remove allBedgraphs folder, comment out if needed
rm -r "$SAVEPATH"allBedgraphs/

echo
echo "Reads copied for $REGIONFILE"
