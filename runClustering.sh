#!/bin/bash

#INPUT: The directory with input file (extension *.out*) = $1
#OUTPUT: The directory where the clustering file is written = $2 (Two files are written: 1. With clustering stats 2. With the read sequence data for each cluster)

#The clustering requires scikit-learn (uses its GMM implementation for clustering)
module load scikit-learn

t=`pwd`
prefix="pacRun."
ext=".out"
ext2=".out.period.sort"

#TO cat the files from the main run directory:
cat $1/*out* > $1/$prefix$ext #the input directory is supplied by the user
#Remove all the lines with "query" and "inf" tags:
grep -v "query" $1/$prefix$ext | grep -v "inf" | sort -nk14 -nk15 -nk16 | awk 'NF' > $1/$prefix$ext2

#Once filtered, the pacMonStr output file can then be used for clustering:

#PARAMETERS for running clustering:
cSep="2.0" #The separation between the means of clusters
probB="0.001" #The binomial probability described by the distribution of reads into two clusters
minReads="4" #The minimum number of reads required to span a TR event for clustering
maxReads="35" #The maximum number of reads required to span a TR event for clustering. Number more than that is probabilistically zero based on data coverage 
outD="clusterRun" #Directory to write the clustering output file
mkdir $2/$outD #The output directory supplied by the user

python $t/gmmCluster/clusteringGMM_AIC_cSep.py $1/$prefix$ext2 $2/$outD 0.001 4 35 $cSep
