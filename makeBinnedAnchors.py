!/usr/bin/env python

import sys
"""
The objective of this tiny script is to generate the binned
anchors file from a blasr m5 alignment and the -ngs output
from tandem repeat finder.

"""

print "format is ./python/makeBinnedAnchors.py reads.m5 trf.reads.ngs"
m5In = sys.argv[1]
trfIn = sys.argv[2]

m5Dict = {}
trfDict = {}

for entry in open(m5In):
    l = entry.strip().split()
    m5Dict[l[0]] = l

for entry in open(trfIn):
    l = entry.strip()
    if l[0] == "@":
        read = l[1:]
        continue
    l = entry.strip().split()
    #repeat events are indexed by their start in
    #    their respective read
    trfDict.setdefault(read,[]).append(l)

for name, entries in trfDict.iteritems():
    m = m5Dict[name]
    for e in entries:
        newEntry = ' '.join(map(str,[m[0],m[2],m[3],int(m[2])+int(e[0]),int(m[2])+int(e[1]),e[2],e[3],e[13],e[0],int(e[1])-int(m[2])]))
        print newEntry

        
