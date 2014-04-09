#!/usr/bin/env python

#####################################################################################################################################################################
#REQUIRES: scikit-learn/scipy
#sys.argv[1] = Sorted file with TR events 
#sys.argv[2] = output directory location 
#python gmmBicClusteringIter.py <INPUT_FILE> <OUTPUT_DIR> 
#####################################################################################################################################################################
import sys
import os
from itertools import *
import numpy as np
from gmmCluster2AlleleAIC_BIC_AD import *

class outL:
    def __init__(self):
        self.Query_id, self.target_id, self.target_start, self.target_end = None,None,0,0
        self.tRSeq, self.dpNum, self.RefNum, self.naiveTR, self.numTR = None, 0.0, 0.0, 0.0, 0.0
        self.ReftRStart, self.ReftREnd, self.RefPeriod = 0,0,0
        self.samples = [] 
        self.Seq = None
        self.tRinQSeq = None
        self.preFlank, self.sufFlank = None, None
        self.ScrTrIn, self.ScrTrOut = 0, 0
        self.preAcc,self.sufAcc,self.trAcc,self.tR_diff_Ex_S = 0.0,0.0,0.0,0

def parseOutFile(file):
    for line in open(file):
        values = line.rstrip("\n").split("[(")[0].split(" [")[1].split("] ")[0]
        #print values
        sampleList = [x for x in values.split(",")]
        values2 = line.rstrip("\n").split(" [")[0].split()
        #print values2
        hit = outL()
        hit.Query_id, hit.target_id, hit.target_start, hit.target_end = values2[0],values2[1],int(values2[2]),int(values2[3])
        hit.preAcc,hit.sufAcc = float(values2[7]), float(values2[8])
        hit.trAcc, hit.tR_diff_Ex_S = float(values2[9]), int(values2[10])
        hit.tRSeq, hit.dpNum, hit.RefNum, hit.naiveTR, hit.numTR = values2[4], float(values2[5]), float(values2[6]), float(values2[16]),float(values2[17])
        hit.ReftRStart, hit.ReftREnd, hit.RefPeriod = int(values2[13]),int(values2[14]),float(values2[15])
        hit.tRinQSeq = ''.join([x if x != '-' else '' for x in values2[12]])
        hit.samples = sampleList

        values3 = line.rstrip("\n").split()[-2:]
        hit.preFlank, hit.sufFlank = values3[0], values3[1]
        hit.Seq = hit.preFlank + hit.tRinQSeq + hit.sufFlank

        maxScrPerQBase = 1.0*(5.0*hit.RefPeriod*hit.dpNum)
        #tR_ScoreDiff = hit.ScrTrOut - hit.ScrTrIn
        alignScrPerQBase = -1.0*(hit.tR_diff_Ex_S) + 5.0 #as ScrTrOut < ScrTrIn. 5.0 is added as number of bases would be 1 more than the difference of tr scores in and out.
        ScoreRatio = alignScrPerQBase/float(maxScrPerQBase)
        #print ScoreRatio
        if (ScoreRatio >= 0.45): #accuracy more than 75% relates to this ratio being greater than 0.5 (check notes for derivations)
           yield hit
        #yield hit
        #if len(hit.samples) == 50:
        #   yield hit

def getReads(hits,statSamples,outfile):
    weight = []
    resultStr = ''
    weight0 = 0
    weight1 = 0
    for j in range(0,len(hits)):
        label0 = 0
        label1 = 0
        temp = []
        for k in range(0,len(statSamples)):
            if statSamples[k][5][j] == 0:
               label0 += 1
            else:
               label1 += 1
        if label0 > label1:
           temp = [hits[j].ReftRStart,hits[j].ReftREnd,hits[j].Query_id,0,hits[j].tRSeq,hits[j].Seq,hits[j].dpNum, hits[j].RefNum, hits[j].naiveTR, hits[j].numTR]
           resultStr = resultStr + "\n" + ' '.join(map(str,temp)) 
           weight0 += 1
        else:
           temp = [hits[j].ReftRStart,hits[j].ReftREnd,hits[j].Query_id,1,hits[j].tRSeq,hits[j].Seq,hits[j].dpNum, hits[j].RefNum, hits[j].naiveTR, hits[j].numTR]
           resultStr = resultStr + "\n" + ' '.join(map(str,temp)) 
           weight1 += 1
    #writeToFile(resultStr,outfile)
    weight = [weight0/float(weight0+weight1),weight1/float(weight0+weight1)]
    return weight 

def writeReadsToFile(hits,call,label,outfile):
    resultStr = ''
    if call == 1:
       for j in range(0,len(hits)):
           temp = [hits[j].ReftRStart,hits[j].ReftREnd,hits[j].Query_id,hits[j].tRSeq,hits[j].tRinQSeq,hits[j].RefNum,hits[j].numTR,hits[j].Seq,0]
           resultStr = resultStr + "\n" + ' '.join(map(str,temp))
    if call == 2:
       for j in range(0,len(hits)):
           temp = [hits[j].ReftRStart,hits[j].ReftREnd,hits[j].Query_id,hits[j].tRSeq,hits[j].tRinQSeq,hits[j].RefNum,hits[j].numTR,hits[j].Seq,label[j]]
           resultStr = resultStr + "\n" + ' '.join(map(str,temp))
    writeToFile(resultStr,outfile)
            
def getProbOfAllele(hits,outfile,c,coV):
    lenSampleList = len(hits[0].samples)
    clusterSamples,statSamples,numAlleleCalls,simMats = [],[],[],[]
    mean1,covar1,mean21,covar21,mean22,covar22 = 0.0,0.0,0.0,0.0,0.0,0.0
    for j in range(0,lenSampleList):
        temp = []
        for k in range(0,len(hits)):
            temp.append(hits[k].samples[j])
        clusterSamples.append(temp)
    for sample in clusterSamples:
        oneAllele,weight,labelD,sM = cluster_TRsBIC(sample,c,coV)
        #print "oneAllele: ",oneAllele, "\nweight:\n", weight,"\nlabelD:\n",labelD,"\nsM:\n",sM
        clustResult = [oneAllele,weight,labelD,sM]
        numAlleleCalls.append(oneAllele)
        statSamples.append(clustResult)
        if sM.shape != (1,):
           simMats.append(sM) #if not an empty simi matrix
    hetz, homz = 0, 0
    for flag in numAlleleCalls:
        if flag == False:
           hetz += 1
        else:
           homz += 1
    if (hetz == 0) and (homz == 0):
       probHomz, probHetz = 0.0, 0.0
       weight = [0.0,0.0]
       numAllele = 0
       probAllele = 0.0
       probBN = 0.0
    else:
       probHomz = homz/float(homz+hetz)
       probHetz = hetz/float(homz+hetz)
       if probHomz >= probHetz:
          numAllele = 1
          weight = [1.0,0.0]
          probAllele = probHomz
          probBN = 1.0
       if probHomz < probHetz:
          labelsHits,LenCl0,LenCl1 = getLabelsFromSimMatrix(simMats,hits)
          #print "labelsHits: \n",labelsHits,"\nLenCl0: ", LenCl0," LenCl1: ", LenCl1
          netNum = LenCl0+LenCl1
          weight = [LenCl0/float(netNum),LenCl1/float(netNum)]
          numAllele = 2
          probAllele = probHetz
          numCL = [weight[0]*netNum,weight[1]*netNum]
          probBN = getProbFromBinom(numCL)
    return numAllele,probAllele,weight,probBN

def getLabelsFromSimMatrix(sML,hits):
    labelSM = []
    labelHits = []
    lenSML = len(sML)
    #dimI, dimJ = sML[0].shape[0], sML[0].shape[1]
    dimJ = sML[0].shape[0]
    #consLabel = np.zeros((dimI,dimJ),dtype=np.int)
    consLabel = np.zeros((dimJ),dtype=np.int)
    #for i in range(0,dimI):
    #for j in range(i+1,dimJ):
    for j in range(1,dimJ):
        for k in range(0,lenSML):
            if sML[k][j] == 1:
               consLabel[j] += 1   
            if sML[k][j] == -1:
               consLabel[j] += -1   
            #if sML[k][i][j] == 1:
            #   consLabel[i][j] += 1   
            #if sML[k][i][j] == -1:
            #   consLabel[i][j] += -1   
    cluster0 = []
    cluster1 = []
    clusterReject = []
    cluster0.append(0)
    #print "consLabel:\n",consLabel
    for j in range(1, dimJ):     
        if consLabel[j] > int(lenSML/2.0):
           cluster0.append(j)
        if consLabel[j] >= 0 and consLabel[j] <= int(lenSML/2.0):
           cluster1.append(j)
        if consLabel[j] < 0:
           clusterReject.append(j)
    #print "cluster0: ",cluster0, " cluster1: ",cluster1, " clusterReject: ",clusterReject
    lenCl0 = len(cluster0)
    lenCl1 = len(cluster1)
    lenClR = len(clusterReject)
    for i in range(0,len(hits)):
        if i in cluster0:
           labelHits.append(0)
        if i in cluster1:
           labelHits.append(1)
        if i in clusterReject:
           labelHits.append(-1)
    return labelHits,lenCl0,lenCl1 

def getClusterPerMethod(naiveTRs,c,hits,coV,outfile2):
    numTRs = len(naiveTRs)
    if numTRs >= 4: #Atleast 4 reads required
       oneAllele,weight,labelD,sM,A2Stats,lb,statsCall = cluster_TRs(naiveTRs,c,coV)
       #oneAllele,weight,labelD,sM,A2Stats,lb = cluster_TRsBIC(naiveTRs,c,coV)
       if oneAllele == False:
          call = 2 #"Heteroz"
          simMats = [sM]
          labelsHits,LenCl0,LenCl1 = getLabelsFromSimMatrix(simMats,hits)
          #print "labelsHits: \n",labelsHits,"\nLenCl0: ", LenCl0," LenCl1: ", LenCl1
          netNum = LenCl0+LenCl1
          weight = [LenCl0/float(netNum),LenCl1/float(netNum)]
          numCL = [weight[0]*netNum,weight[1]*netNum]
          probBin = getProbFromBinom(numCL)
       else:
          call = 1 #"Homoz"
          weight = [1.0,0.0]
          probBin = 1.0
          labelsHits = []
       writeReadsToFile(hits,call,labelsHits,outfile2)
       return call,weight,probBin,A2Stats,lb,statsCall

def writeToFile(resultStr,outfile):
    outfile.write("%s\n" %(resultStr))
    outfile.flush()

def binQtrs():
    inm5fn = sys.argv[1] #input .out file with  sorted sample data
    tempm5 = inm5fn.split("/")[-1].split(".")[0]
    outBinfn = sys.argv[2] #out directory with the clustering results: Probability of cluster1 or cluster2
    probCutOff = float(sys.argv[3]) #A PARAMETER for BN calculations
    minNumReads = int(sys.argv[4])
    maxNumReads = int(sys.argv[5])
    coV = 0.0175 #A PARAMETER for minCovariance value for GMM
    cSep = float(sys.argv[6]) #A PARAMETER for checking the c-separation between the means of two clusters

    #fnExt = str(coV)
    fnExt2 = str(cSep)
    hgT_outFn = outBinfn + "/" + tempm5 + "_" + '_'.join(map(str,str(probCutOff).split("."))) + "_" + str(minNumReads) + "_minCovar_0175_" + fnExt2 + ".clusterCalls_TR.txt"
    outfile = open(hgT_outFn, 'w')
    #outfile.write("tRstart tRend refPeriod probCutOff minNumReads TotalReads allele_V probBN_V allele_D probBN_D allele_H probBN_H\n")
    outfile.write("tRstart tRend refPeriod refNum coV probCutOff minNumReads TotalReads allele_H cl1mean cl2mean probBN_H ADpH cSepH\n")
    hgT_outFn2 = outBinfn + "/" + tempm5 + "_" + '_'.join(map(str,str(probCutOff).split("."))) + "_" + str(minNumReads) + "_clusterReadData_TR.txt"
    outfile2 = open(hgT_outFn2, 'w')

    hits = parseOutFile(inm5fn)
    for k,g in groupby(hits, key=lambda x: (x.ReftRStart,x.ReftREnd)):
      for k1,g1 in groupby(g, key=lambda x: x.RefPeriod):
        hits = list(g1)
        localDict = {}
        for hit in hits:
            localDict[hit.Query_id] = hit
            #print hit.Query_id
        #print "length of localDict", len(localDict)
        #qidList = localDict.keys()
        #print qidList
        if len(localDict) >= minNumReads and len(localDict) <= maxNumReads:
           pHmmTRs = []
           for qid,outObj in localDict.iteritems():
               numTRlocal = outObj.numTR
               pHmmTRs.append(numTRlocal)
           #callpH,weightpH,probBinpH,ADpH,cSepH,statsCall = getClusterPerMethod(pHmmTRs,cSep,hits,coV,outfile2)
           localList = localDict.values()
           callpH,weightpH,probBinpH,ADpH,cSepH,statsCall = getClusterPerMethod(pHmmTRs,cSep,localList,coV,outfile2)
           statA = np.around(np.asarray(statsCall),decimals=2)
           resultList = [hits[0].ReftRStart,hits[0].ReftREnd,hits[0].RefPeriod,hits[0].RefNum,coV,probCutOff,minNumReads,len(localDict),callpH,statA[0],statA[1],probBinpH,ADpH,cSepH]
           resultStr = ' '.join(map(str,resultList))
        else:
           resultStr = str(hits[0].ReftRStart) + " " + str(hits[0].ReftREnd) + " Less than 4 or more than 35 reads spanning TR event"
        writeToFile(resultStr,outfile)
    outfile.close() 
    outfile2.close() 

binQtrs()
