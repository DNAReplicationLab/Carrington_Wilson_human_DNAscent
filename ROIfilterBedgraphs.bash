#!/bin/bash

# Takes detect file and bedgraphs from dnascent and copies bedgraphs corresponding to region of interest to a new folder

# usage bash ROIfilterBedgraphs.bash -d </path/to/.detect/file> -b </path/to/bedgraphs> -s </save/path> -c <chrom> -a <start> -e <end> -n <ROI_name> optional [-w <value of ROI extension> ]]
# -c <chrom> -a <start> -e <end> -n <origin_name> are all for your region of interest
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
		-c )	shift
			CHROM="$1" ;;
		-a )	shift
			START="$1" ;;
		-e )	shift
			END="$1" ;;
		-n )	shift
			NAME="$1" ;;
		-w )	shift
			WIDE="$1" ;;
	esac
		shift
done

# extend region if -w
if [ "$WIDE" != 0 ]; then
	echo
	echo "Region will be extended by $WIDE each side."
	START=$(expr "$START" - "$WIDE")
	END=$(expr "$END" + "$WIDE")
fi

# print variables
echo
echo "Reads to be copied to $SAVEPATH for $NAME, $CHROM, $START, $END from $BEDGRAPHS using $DETECT"

#make list of reads that fall within ROI + 100kb either side
mkdir "$SAVEPATH"temp
grep "^>" "$DETECT" > "$SAVEPATH"temp/readNamesDetect.txt
cat "$SAVEPATH"temp/readNamesDetect.txt | grep "$CHROM" > "$SAVEPATH"temp/readNamesChr.txt
#pass start and end coordinates to awk, keep lines where reads are within ROI genome coordinates.
awk -v start="$START" -v end="$END" 'start < $4 && $3 < end' "$SAVEPATH"temp/readNamesChr.txt > "$SAVEPATH""$NAME".readNames.txt
cat "$SAVEPATH""$NAME".readNames.txt | cut -f2 -d ">" | cut -f1 -d " " > "$SAVEPATH"temp/readNamesROI.txt
echo
echo "Reads identified for $CHROM, start $START, end $END"

#copy all bedgraphs to single folder
mkdir "$SAVEPATH"temp/allBedgraphs "$SAVEPATH""$NAME".bedgraphs
cp "$BEDGRAPHS"*/*.bedgraph "$SAVEPATH"temp/allBedgraphs/
cd "$SAVEPATH"temp/allBedgraphs/

#copy ROI bedgraphs to ROI folder
cp $(ls . | grep -f "$SAVEPATH"temp/readNamesROI.txt) "$SAVEPATH""$NAME".bedgraphs/
cd "$WD"

#remove temp folder, comment out line below if want to check intermediate files
rm -r "$SAVEPATH"temp/

echo
echo "Reads copied for $NAME, $CHROM, start $START, end $END"
