import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

target = 'III'
length = 2452883
BrdUCalls = [0]*length
coverage = [0]*length
BrdUCalls = np.array(BrdUCalls)
coverage = np.array(coverage)

#import detect data
f = open(sys.argv[1],'r')
for line in f:

        if line[0] == '#':
    
            continue
        
        
        if line[0] == '>':

		splitLine = line.rstrip().split(' ')
		chromosome = splitLine[1]

		continue

	else:

		splitLine = line.rstrip().split('\t')

		if chromosome == target:

			coverage[int(splitLine[0])] += 1
			if float(splitLine[1]) > 0.12:

				BrdUCalls[int(splitLine[0])] += 1



xBrdU = []
yBrdU = []
for i in range( 0, length, 1000 ):

	if float(sum( coverage[i:i+1000])) == 0.0:
		continue
	else:
		yBrdU.append(float(sum( BrdUCalls[i:i+1000] )) / float(sum( coverage[i:i+1000])))
		xBrdU.append(i+500)

yBrdUSmooth = np.convolve(yBrdU, np.ones((10,))/10, mode='same')


f = open('BrdU_data_chrIII.plot','w')
f.write(str('variableStep chrom=chr3' + '\n'))
for i, v in enumerate(yBrdUSmooth):
    f.write(str(xBrdU[i]) + '\t' + str(v) + '\n')
f.close()
