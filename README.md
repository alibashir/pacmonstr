PacmonsTR
=========

#### Tandem Repeat Detection for Long Read Sequence Data

The pipeline works in the following steps:

1. The initial run takes in an alignment file (alignment of long reads with the reference, such as, hg19) and a reference Tandem Repeat table as INPUT and generates an OUTPUT file with estimated tanderm repeat region in the Long read, expected value of the number of tandem repeat elements, estimated Structural variation region in the tandem repeat.

2. The output file generated from Step A, is used (after processing) for clustering. The clustering uses GMM based clustering method, which uses AIC as the model selection criterion, coupled with c-separation. Allelic call (homozygous/heterozyous) is made based on the clustering criteria and it outputs two files. One file has all the calls and the other has the read data for each cluster.

3. The read data for each cluster is taken and a consensus sequence is generated. 

4. The consensus sequence can be taken as an input to StepA and re-estimation of TR multiplicity can be performed.
