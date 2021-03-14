#----------------------------------------------------------
# Copyright 2020-2021 Earlham Institute
# Written by Conrad Nieduszynski (conrad.nieduszynski@earlham.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

# Usage: python DetectToEnsembleWig.py path/to/input/detect/file window_in_kb path/for/output/wig/file
# example command:
# python DetectToEnsembleWig.py test.detect 100 test.wig

import sys
import numpy as np

chromosomes = 24
maxLength = 248956422 // 1000

BrdUCalls = np.zeros(shape=(chromosomes,maxLength))
coverage = np.zeros(shape=(chromosomes,maxLength))
threshold = 0.5
SmallWindow = 1000
LargeWindow = int(sys.argv[2])

#import detect data
f = open(sys.argv[1],'r')
for line in f:

	if line[0] == '#':

		continue
	
	
	if line[0] == '>':

		splitLine = line.rstrip().split(' ')
		chromosome = splitLine[1]
# 		print(chromosome)
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
		if float(splitLine[1]) > threshold:

			BrdUCalls[MyChromosome-1][int(splitLine[0])//SmallWindow] += 1

#export the ensemble data as a wig file
ensembleBrdU = np.zeros(shape=(chromosomes,maxLength))

f = open(sys.argv[3],'w')
for chromo in range(len(ensembleBrdU)):
	if chromo < 23:
		f.write('variableStep chrom=chr{}'.format(chromo+1) + '\n')
	elif chromo == 23:
		f.write('variableStep chrom=chrX' + '\n')
	elif chromo == 24:
		f.write('variableStep chrom=chrY' + '\n')
	
	for i in range( 0, maxLength, LargeWindow ):

		if float(sum( coverage[chromo][i:i+LargeWindow])) == 0.0:
			continue
		else:
			ensembleBrdU[chromo][i + LargeWindow//2] = float(sum( BrdUCalls[chromo][i:i+1000] )) / float(sum( coverage[chromo][i:i+1000]))
			f.write(str((i + LargeWindow//2)*1000) + '\t' + str("%.3f" % ensembleBrdU[chromo][i + LargeWindow//2]) + '\n')

f.close()

# yBrdUSmooth = np.convolve(yBrdU, np.ones((10,))/10, mode='same')
# 
# 
# f = open('BrdU_data_chrIII.plot','w')
# f.write(str('variableStep chrom=chr3' + '\n'))
# for i, v in enumerate(yBrdUSmooth):
#     f.write(str(xBrdU[i]) + '\t' + str(v) + '\n')
# f.close()
