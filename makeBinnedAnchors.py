#!/usr/bin/env python

import sys
import bisect
"""
The objective of this tiny script is to generate the binned
anchors file from a blasr m5 alignment and a tandem repeat reference file.
The presumed fields in the tandem repeat file are 

593	chr1	1108151	1108213	trf	23	2.7	22	79	13	63	75	1	0	22	0.88	AAAATAAAATATAAAAAAAATA
"""

print "format is ./python/makeBinnedAnchors.py reads.m5 trf.reads.ngs"
m5In = sys.argv[1]
repeatFile = sys.argv[2]
flank= int(sys.argv[3])

def searchTRs(TRlist, start, stop, flank=100):
    if stop - start < 2* flank:
        #read too short to contain TR and flank sequences
        return -1
    flankLeft = start + flank
    flankRight = stop - flank
    i1 =bisect.bisect(TRlist,(flankLeft,flankLeft+1))
    candidates1 = TRlist[i1:]
    i2 = bisect.bisect(candidates1,(flankRight,flankRight))
    return candidates1[:i2] #a list of all TRs in a given read s.t. each TR is flanked by at least `flank` bp of seq

m5Dict = {}
repeatDict = {}
repeatsFull = {}    
for entry in open(repeatFile):
    l = entry.strip().split()
    ch,start,stop = l[1:4]
    start = int(start)
    stop = int(stop)
    repeatDict.setdefault(ch,[]).append((start,stop)) #this one is for binomially searching
    repeatsFull.setdefault(ch,{})[(start,stop)] = l #this one is for retrieving the full entry

#sort each chromosome by start position of tuples
# these can then be searched using bisection to find an overlapping TR given the (start,stop) of a given read
print 'Sort chromosomes by TR start position, then TR end position.'

for key in repeatDict.keys():
    repeatDict[key].sort()
    print key,'sorting complete.'

print 'Done sorting.'

for entry in open(m5In):
    fOut = open(m5In.replace('.m5','.binned_anchors'),'w')
    m = entry.split()
    if m[5] not in repeatDict:
	continue
    TRsInEntry = searchTRs(repeatDict[m[5]], int(m[7]),int(m[8]),flank)
    if TRsInEntry == -1:
	print m
        continue
    if len(TRsInEntry)==0:
	continue
    for tr in TRsInEntry:
        e=repeatsFull[m[5]][tr]
        fOut. write(' '.join(map(str,[m[0],m[7],m[8],e[2], e[3],e[5],"{:.2f}".format(float(e[15])),e[16],int(e[2]) -int(m[7]),int(m[8]) -int(e[3])]))+'\n')
        print ' '.join(map(str,[m[0],m[7],m[8],e[2], e[3],e[5],"{:.2f}".format(float(e[15])),e[16],int(e[2]) -int(m[7]),int(m[8]) -int(e[3])]))

