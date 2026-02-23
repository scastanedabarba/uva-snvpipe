#!/usr/bin/env python

import sys, os
import argparse

### accept arguments
parser = argparse.ArgumentParser("")
parser.add_argument("-t", type=str, help="MASH-triangle-results", required=True)
parser.add_argument("-o", type=str, help="MASH-closest-ref", required=True)
args = parser.parse_args()
args = vars(args)

## first read the matrix to get list of references and all sequence IDs
myphyliporder = []
myrefs = []
mysamples = []
fh = open(args['t'])
fh.readline()
for line in fh:
	line  = line.strip()
	arr = line.split("\t")
	arr = [x.strip() for x in arr]

	myphyliporder += [arr[0].split("/")[-1]]	
	
	## list of reference accessions
	if "GCF_" in arr[0]:
		myref = arr[0].split("/")[-1]
		if myref not in myrefs:
			myrefs += [myref]
	else:
		mysample = arr[0].split("/")[-1]
		if mysample not in mysamples:
			mysamples += [mysample]

fh.close()


## read the matrix again to get the pairwise distances
pairwisedist = dict()
fh = open(args['t'])	
fh.readline()
for line in fh:
	line  = line.strip()
	arr = line.split("\t")
	arr = [x.strip() for x in arr]	
	
	for i in range(1,len(arr)):
		k1 = arr[0].split("/")[-1]

		#k2 = mylinelist[i-1]	
		k2 = myphyliporder[i-1]
		pairwisedist.setdefault(k1,{}).setdefault(k2,0.0)
		pairwisedist[k1][k2] = float(arr[i])
fh.close()


### calculate the overall distance of each reference for all samples
refdist = dict()
for ref in myrefs:
	refdist.setdefault(ref, 0.0)
	for sample in mysamples:
		try:
			refdist[ref] += pairwisedist[ref][sample]
		except:
			refdist[ref] += pairwisedist[sample][ref]


### get reference with shortest overall distance to all isolates
sorted_mydistances = dict(sorted(refdist.items(), key=lambda x:x[1], reverse=False))
myref = list(sorted_mydistances.items())[0]

### Write output linelist
ofh = open(args['o'], "w")
myref = [str(x) for x in myref]
print("\t".join(myref), file=ofh)
ofh.close()
sys.exit()
