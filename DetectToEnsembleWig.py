#----------------------------------------------------------
# Copyright 2020-2021 Earlham Institute
# Written by Conrad Nieduszynski (conrad.nieduszynski@earlham.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

# Usage: python DetectToEnsembleWig.py path/to/input/detect/file window_in_kb path/for/output/files path/to/nascent/bed/file
# example command:
# python DetectToEnsembleWig.py test.detect 100 output_file_directory nascent.bed
# requires python 3.3+

import os
import sys
import numpy as np

chromosomes = 24
maxLength = 248956422 // 1000

BrdUCalls = np.zeros(shape=(chromosomes,maxLength))
nascentBrdUCalls = np.zeros(shape=(chromosomes,maxLength))
coverage = np.zeros(shape=(chromosomes,maxLength))
nascentCoverage = np.zeros(shape=(chromosomes,maxLength))
threshold = 0.5
SmallWindow = 1000
LargeWindow = int(sys.argv[2])
OutputPath = str(sys.argv[3])
BEDfile = str(sys.argv[4])
allReadsBEDfile = os.path.join(OutputPath, "AllReadsEnsemble.bedgraph")
allReadsWIGfile = os.path.join(OutputPath, "AllReadsEnsemble.wig")
nascentReadsBEDfile = os.path.join(OutputPath, "NascentReadsEnsemble.bedgraph")
nascentReadsWIGfile = os.path.join(OutputPath, "NascentReadsEnsemble.wig")
nascentReadIDs = []
nascent = False

print("Will generate file:", allReadsBEDfile, allReadsWIGfile, nascentReadsBEDfile, nascentReadsWIGfile)

try:
    os.makedirs(OutputPath)
except FileExistsError:
    # directory already exists
    pass

print("Reading nascent BED file...")

# open file in read mode
f = open(BEDfile, 'r')
for row in f:
	splitLine = row.rstrip().split('\t')
	try:
		nascentReadIDs.append(splitLine[3])
	except:
		continue
f.close()

print("Reading detect file...")

#import detect data
f = open(sys.argv[1],'r')
for line in f:

	if line[0] == '#':

		continue
	
	
	if line[0] == '>':

		splitLine = line.rstrip().split(' ')
		readID = splitLine[0][1:]
		chromosome = splitLine[1]
		if readID in nascentReadIDs:
			nascent = True
		else:
			nascent = False
		continue


	else:
		
		if chromosome[3:] == 'X':
			MyChromosome = 23
			
		elif chromosome[3:] == 'Y':
			MyChromosome = 24
		
		elif chromosome[3:].isnumeric():
			MyChromosome = int(chromosome[3:])
		
		else:
			continue

		splitLine = line.rstrip().split('\t')
		coverage[MyChromosome-1][int(splitLine[0])//SmallWindow] += 1
		if nascent:
			nascentCoverage[MyChromosome-1][int(splitLine[0])//SmallWindow] += 1
			if float(splitLine[1]) >= threshold:
				nascentBrdUCalls[MyChromosome-1][int(splitLine[0])//SmallWindow] += 1
			
		if float(splitLine[1]) >= threshold:
			BrdUCalls[MyChromosome-1][int(splitLine[0])//SmallWindow] += 1

#export the ensemble data as a wig file
ensembleBrdU = np.zeros(shape=(chromosomes,maxLength))
nascentEnsembleBrdU = np.zeros(shape=(chromosomes,maxLength))

print("Writing files")

f1 = open(allReadsWIGfile,'w')
f2 = open(allReadsBEDfile,'w')
f3 = open(nascentReadsWIGfile,'w')
f4 = open(nascentReadsBEDfile,'w')
for chromo in range(len(ensembleBrdU)):
	if chromo < 22:
		myChromo = chromo + 1
		f1.write('variableStep chrom=chr{}'.format(chromo+1) + '\n')
	elif chromo == 22:
		myChromo = 'X'
		f1.write('variableStep chrom=chrX' + '\n')
	elif chromo == 23:
		myChromo = 'Y'
		f1.write('variableStep chrom=chrY' + '\n')
	
	for i in range( 0, maxLength, LargeWindow ):

		if float(sum( coverage[chromo][i:i+LargeWindow-1])) == 0.0:
			continue
		else:
			ensembleBrdU[chromo][i + LargeWindow//2] = float(sum( BrdUCalls[chromo][i:i+LargeWindow-1] )) / float(sum( coverage[chromo][i:i+LargeWindow-1]))
			f1.write(str((i + LargeWindow//2)*1000) + '\t' + str("%.3f" % ensembleBrdU[chromo][i + LargeWindow//2]) + '\n')
			f2.write('chr' + str(myChromo) + '\t' + str((i - LargeWindow//2)*1000) + '\t' + str((i + LargeWindow//2)*1000) + '\t' + str("%.3f" % ensembleBrdU[chromo][i + LargeWindow//2]) + '\n')
f1.close()
f2.close()

f3 = open(nascentReadsWIGfile,'w')
f4 = open(nascentReadsBEDfile,'w')
for chromo in range(len(nascentEnsembleBrdU)):
	if chromo < 22:
		myChromo = chromo + 1
		f3.write('variableStep chrom=chr{}'.format(chromo+1) + '\n')
	elif chromo == 22:
		myChromo = 'X'
		f3.write('variableStep chrom=chrX' + '\n')
	elif chromo == 23:
		myChromo = 'Y'
		f3.write('variableStep chrom=chrY' + '\n')
	
	for i in range( 0, maxLength, LargeWindow ):

		if float(sum( nascentCoverage[chromo][i:i+LargeWindow-1])) == 0.0:
			continue
		else:
			nascentEnsembleBrdU[chromo][i + LargeWindow//2] = float(sum( nascentBrdUCalls[chromo][i:i+LargeWindow-1] )) / float(sum( nascentCoverage[chromo][i:i+LargeWindow-1]))
			f3.write(str((i + LargeWindow//2)*1000) + '\t' + str("%.3f" % nascentEnsembleBrdU[chromo][i + LargeWindow//2]) + '\n')
			f4.write('chr' + str(myChromo) + '\t' + str((i - LargeWindow//2)*1000) + '\t' + str((i + LargeWindow//2)*1000) + '\t' + str("%.3f" % nascentEnsembleBrdU[chromo][i + LargeWindow//2]) + '\n')

f3.close()
f4.close()


# yBrdUSmooth = np.convolve(yBrdU, np.ones((10,))/10, mode='same')
# 
# 
# f = open('BrdU_data_chrIII.plot','w')
# f.write(str('variableStep chrom=chr3' + '\n'))
# for i, v in enumerate(yBrdUSmooth):
#     f.write(str(xBrdU[i]) + '\t' + str(v) + '\n')
# f.close()
