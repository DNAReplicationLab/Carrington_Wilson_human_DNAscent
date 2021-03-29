#----------------------------------------------------------
# Copyright 2020-2021 Earlham Institute
# Written by Conrad Nieduszynski (conrad.nieduszynski@earlham.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

# Usage: python DetectSummary.py path.to.detect.file 
# example command:
# python DetectSummary.py test_output.detect 

import sys
import numpy as np
from scipy.ndimage.filters import uniform_filter1d
from scipy.signal import find_peaks
import copy

detectThreshold = 0.5
readThreshold = 0.05
first = True
Window = 290
sWindow = Window * 3
lWindow = sWindow * 5
stepBuff = []
maxStepBuff = []
seStepBuff = []
readCount = 0
nascentCount = 0
nascBEDbuff = []
parBEDbuff = []

headerText = ['chr',
   'start',
   'end',
   'read',
   '[BrdU]',
   'strand',
   'mapped length',
   'highest BBT',
   'largest step',
   'start-end diff']

def processRead():
	global nascentCount
	BrdU_array = np.array(detectBuff)
	coorArray = np.array(coorBuff)
	sWindowDiff = uniform_filter1d(BrdU_array, size = sWindow, mode = 'wrap', origin = (sWindow//2)-1) - uniform_filter1d(BrdU_array, size = sWindow, mode = 'wrap', origin = -((sWindow//2)-1))
	absWindowDiff = abs(sWindowDiff)
	readBBT = np.sum(BrdU_array >= detectThreshold)/BrdU_array.size
	hardB = np.floor(BrdU_array / detectThreshold)
	BBT = uniform_filter1d(hardB, size = Window, origin = 0)
	if readBBT >= readThreshold and BrdU_array.size > (2*sWindow):
		stepIndices, _ = find_peaks(absWindowDiff[sWindow:-sWindow], distance = lWindow)
		stepCoor = coorArray[stepIndices + sWindow]
		stepSize = absWindowDiff[stepIndices + sWindow]
		nascentCount += 1
		for step in sWindowDiff[sWindow:-sWindow]:
			stepBuff.append(float(step))
		endStep = absWindowDiff[len(absWindowDiff) - 1]
		if stepSize.size > 0:
		   maxStepBuff.append(max(stepSize))
		else:
		   maxStepBuff.append(0)
		seStepBuff.append(endStep)
		BEDtemp = (chromosome,
		   mappingStart,
		   mappingEnd,
		   readID,
		   str("%.3f" % (readBBT)),
		   strand,
		   (mappingEnd - mappingStart),
		   str("%.3f" % (BBT[Window:-Window].max())),
		   str("%.3f" % (absWindowDiff[sWindow:-sWindow]).max()),
		   str("%.3f" % (endStep)),
		   stepCoor, stepSize)
		nascBEDbuff.append(BEDtemp)
	elif BrdU_array.size > (2*sWindow):
		BEDtemp = (chromosome,
		   mappingStart,
		   mappingEnd,
		   readID,
		   str("%.3f" % (readBBT)),
		   strand,
		   (mappingEnd - mappingStart),
		   str("%.3f" % (BBT[Window:-Window].max())),
		   str("%.3f" % (abs(sWindowDiff[sWindow:-sWindow]).max())),
		   str("%.3f" % (abs(sWindowDiff[len(sWindowDiff) - 1]))),)
		parBEDbuff.append(BEDtemp)


#import detect data
f = open(sys.argv[1],'r')
for line in f:

	if line[0] == '#':

		continue
	
	
	if line[0] == '>':

		readCount += 1

		if not first:
			processRead()

		splitLine = line.rstrip().split(' ')
		readID = splitLine[0][1:]
		chromosome = splitLine[1]
		mappingStart = int(splitLine[2])
		mappingEnd = int(splitLine[3])
		readStrand = splitLine[4]
		if readStrand == 'fwd':
			strand = '+'
		else:
			strand = '-'
				
		coorBuff = []
		detectBuff = []
		first = False

		continue

	else:

		splitLine = line.rstrip().split()
#		print(splitLine[0], splitLine[1])
		coorBuff.append(int(splitLine[0]))
		detectBuff.append(float(splitLine[1]))

#now process the last read
processRead()

stepArray = np.array(stepBuff)
stepSD = stepArray.std()
stepMean = stepArray.mean()
maxStepArray = np.array(maxStepBuff)
seStepArray = np.array(seStepBuff)
nascBEDarray = np.array(nascBEDbuff, dtype=object)
nascBEDarray[:, 11] = nascBEDarray[:, 11] / stepSD
nascStepArray = np.concatenate(nascBEDarray[:, 11], axis=None)
parBEDarray = np.array(parBEDbuff)
print(nascStepArray[:10])

print('## Total reads:', readCount)
print('## Nascent reads:', nascentCount, ' ( Total BrdU fraction in read >', readThreshold, ')')
print('## Total length of nascent reads:', nascBEDarray[:, 6].sum(axis=0))
print('##')
print('## STD in steps:', str("%.3f" % (stepArray.std())), ' (', sWindow, 'T positions )')
print('## Mean in steps:', str("%.3f" % (stepArray.mean())), ' (', sWindow, 'T positions )')
print('## Max step size:', str("%.3f" % (abs(stepArray).max())), ' (', sWindow, 'T positions )')
print('## Min step size:', str("%.3f" % (abs(stepArray).min())), ' (', sWindow, 'T positions )')
print('##')
#These calculations are relative to zero, but they should be relative to the mean
#The mean is very close to zero, so probably won't make much difference
#To correct this I'll need to use the step values before abs
print('## Reads with step > 1 sd:', np.sum(maxStepArray >= stepSD))
print('## Reads with step > 2 sd:', np.sum(maxStepArray >= 2*stepSD))
print('## Reads with step > 3 sd:', np.sum(maxStepArray >= 3*stepSD))
print('## Reads with step > 4 sd:', np.sum(maxStepArray >= 4*stepSD))
print('##')
print('## Number of steps > 2 sd:', np.sum(nascStepArray >= 2))
print('## Number of steps > 3 sd:', np.sum(nascStepArray >= 3))
print('## Number of steps > 4 sd:', np.sum(nascStepArray >= 4))
print('##')
print('## Reads with start-end diff > 1 sd:', np.sum(seStepArray >= stepSD))
print('## Reads with start-end diff > 2 sd:', np.sum(seStepArray >= 2*stepSD))
print('## Reads with start-end diff > 3 sd:', np.sum(seStepArray >= 3*stepSD))
print('## Reads with start-end diff > 4 sd:', np.sum(seStepArray >= 4*stepSD))

boolSEstep = seStepArray > 3*stepSD
boolStep = maxStepArray > 3*stepSD
txtHeader = '\t'.join(headerText)

np.savetxt('nascent.txt', nascBEDarray, fmt='%s', header = txtHeader, delimiter='\t')
np.savetxt('parental.txt', parBEDarray, fmt='%s', header = txtHeader, delimiter='\t')

np.savetxt('nascent.bed', nascBEDarray[:, :6], fmt='%s', delimiter='\t')
np.savetxt('parental.bed', parBEDarray[:, :6], fmt='%s', delimiter='\t')

np.savetxt('start_end_diff.txt', nascBEDarray[boolSEstep], fmt='%s', header = txtHeader, delimiter='\t')
np.savetxt('start_end_diff.bed', nascBEDarray[boolSEstep, :6], fmt='%s', delimiter='\t')

np.savetxt('3SD_step.txt', nascBEDarray[boolStep], fmt='%s', header = txtHeader, delimiter='\t')
np.savetxt('3SD_step.bed', nascBEDarray[boolStep, :6], fmt='%s', delimiter='\t')

