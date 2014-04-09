#!/usr/bin/env python

#####################################################################################################################################################################
#####################################################################################################################################################################
import sys
import math
import os
from itertools import *
from sklearn import mixture as M
from scipy.stats import anderson as ad
from scipy.stats import binom as bN
import numpy as np
#import cProfile

def getZfactor(g_stats):
    zF = -1.0
    sTd1 = math.sqrt(g_stats[0][1])
    sTd2 = math.sqrt(g_stats[1][1])
    diffMeans = abs(g_stats[0][0]-g_stats[1][0])
    if diffMeans != 0.0:
       zF = 1.0 - 3.0*(sTd1+sTd2)/diffMeans
    return zF

def getCseparationGauss(g_stats,c):
    lowerBound = []
    st1 = g_stats[0]
    st2 = g_stats[1]
    sTd12 = max(np.sqrt(g_stats[0][1]),np.sqrt(g_stats[1][1]))
    diffMeans = abs(g_stats[0][0]-g_stats[1][0])
    lowerBound.append(diffMeans)
    lowerBound.append(c*sTd12)
    #for I in c:
    #    lowerBound.append(I*sTd12)
    return lowerBound

def getProbFromBinom(numCl):
    minWt = int(min(numCl))
    numData = int(numCl[0] + numCl[1])
    probMinWt = bN.cdf(minWt,numData,0.5)
    return probMinWt

def andersonDarling(X):
    A2stats = ad(X) #anderson-darling EDF test as implemented in scipy, returns unmodified value of the statistics
    numTR = np.size(X)
    #modified value of A2stats
    #modA2 = (1.0 + 4.0/numTR - 25.0/(numTR*numTR))
    modA2 = (1.0 + 0.75/numTR + 2.25/(numTR*numTR)) #for the case where both mean and variance is unknown, D'Agostino
    A2statsMod = A2stats[0]*modA2
    #if A2statsMod < 1.869: #significance level alpha = 0.0001, Ref: G-means
    #if A2statsMod < 1.159: #significance level alpha = 0.005, Ref: D'Agostino (1986) in Table 4.7, p. 123
    #   acceptHypothesis0 = True
    #else:
    #   acceptHypothesis0 = False
    return A2statsMod
 
def removeLabel(TR_list,label,cl,L):
    TR_listR = []
    TR_listLabel = -10
    for i in range(0,len(TR_list)):
        if label[i] != cl:
           TR_listR.append(TR_list[i])
        else:
           TR_listLabel = (i)
    return TR_listR,TR_listLabel

def decouple(TR_listNA,label):
    temp = list(label)
    TR_listNA.reverse()
    for l in TR_listNA:
        temp.insert(l,-1)
    return temp

def getSimMatrix(label):
    lenLab = len(label)
    if lenLab != 0:
       #sM = np.zeros((lenLab,lenLab),dtype=np.int)
       sM = np.zeros((lenLab),dtype=np.int)
       label_i = label[0]
       #for i in range(0,lenLab):
       sM[0] = 1
       for j in range(1,lenLab):
         if label_i != -1 and label[j] != -1:
            if label_i == label[j]:
               sM[j] = 1
            else:
               sM[j] = 0
         else:
            sM[j] = -1
    else:
       sM = np.zeros((1),dtype=np.int)
    return sM

def cluster_TRsBIC(TR_list,stdParam,coV):
    #np.random.seed(1)
    epsilon = 0.01
    flag = True
    TR_listNA = []
    bic = np.zeros((2),dtype=np.float_)
    while (flag):
        #print "TR_list_start of While loop:\n",TR_list
        numTR = len(TR_list)
        TR_Data = np.zeros((numTR,1),dtype=np.float_) #a column matrix
        TR_Data[:,0] = np.array(TR_list)
        if numTR == 2:
           g1_stats = (np.mean(TR_Data),np.var(TR_Data))
           label = np.array([0,0])
           flag = False
           bic[1], bic[0] = 0, 0
           break
        else:
            #Calling GMM routines here...
            g1 = M.GMM(n_components=1,random_state=1,min_covar=coV)
            g1.fit(TR_Data)
            bic[0] = g1.bic(TR_Data)
            g1_stats = (g1.means_[0][0],g1.covars_[0][0])
            g2 = M.GMM(n_components=2,random_state=1,min_covar=coV)
            g2.fit(TR_Data)
            bic[1] = g2.bic(TR_Data)
            label = g2.predict(TR_Data)
            weight = [g2.weights_[0],g2.weights_[1]]
            g2_stats = ((g2.means_[0][0],g2.covars_[0][0]),(g2.means_[1][0],g2.covars_[1][0]))
            if bic[1] < bic[0]:
               if weight[0]*numTR > 1.0+epsilon and weight[1]*numTR > 1.0+epsilon:
                  flag = False
                  break
               else:
                  if weight[0]*numTR <= 1.0+epsilon:
                     TR_list,TR_listLabel = removeLabel(TR_list,label,0,len(TR_listNA)) #remove label 0 and return the TR list
                  elif weight[1]*numTR <= 1.0+epsilon:
                     TR_list,TR_listLabel = removeLabel(TR_list,label,1,len(TR_listNA)) #remove label 1 and return the TR list
                  TR_listNA.append(TR_listLabel)
                  #print "TR_listNA:\n",TR_listNA
            else:
               break
    tempLabelD = decouple(TR_listNA,label)
    #print "g1_stats: ",g1_stats, " g2_stats: ",g2_stats, " bic: ",bic
    #print "tempLabelD :", tempLabelD
    A2stats = andersonDarling(TR_Data[:,0])
    if bic[1] < bic[0]:
       lb = getCseparationGauss(g2_stats,stdParam)
       if lb[0] > lb[1]:
          oneAllele = False
          sM = getSimMatrix(tempLabelD)
          labelD = tempLabelD
          cl0 = [lab for lab in tempLabelD if lab == 0]
          cl1 = [lab for lab in tempLabelD if lab == 1]
          weight = [len(cl0)/float(len(cl0)+len(cl1)),len(cl1)/float(len(cl0)+len(cl1))]
       else:
          oneAllele = True
          sM = np.zeros((1),dtype=np.int)
          weight = [1.0,0.0]
          labelD = [0 for lab in tempLabelD]
    else:
       #print "homoz"
       lb = [0.0,0.0]
       oneAllele = True
       weight = [1.0,0.0]
       sM = np.zeros((1),dtype=np.int)
       labelD = [0 for lab in tempLabelD]
    return oneAllele,weight,labelD,sM,A2stats,lb
 
def cluster_TRs(TR_list,stdParam,coV):
    #np.random.seed(1)
    epsilon = 0.01
    flag = True
    TR_listNA = []
    aic = np.zeros((2),dtype=np.float_)
    while (flag):
        #print "TR_list_start of While loop:\n",TR_list
        numTR = len(TR_list)
        TR_Data = np.zeros((numTR,1),dtype=np.float_) #a column matrix
        TR_Data[:,0] = np.array(TR_list)
        if numTR == 2:
           g1_stats = (np.mean(TR_Data),np.var(TR_Data))
           label = np.array([0,0])
           flag = False
           aic[1], aic[0] = 0, 0
           break
        else:
            #Calling GMM routines here...
            g1 = M.GMM(n_components=1,random_state=1,min_covar=coV)
            g1.fit(TR_Data)
            aic[0] = g1.aic(TR_Data)
            g1_stats = (g1.means_[0][0],g1.covars_[0][0])
            g2 = M.GMM(n_components=2,random_state=1,min_covar=coV)
            g2.fit(TR_Data)
            aic[1] = g2.aic(TR_Data)
            label = g2.predict(TR_Data)
            weight = [g2.weights_[0],g2.weights_[1]]
            g2_stats = ((g2.means_[0][0],g2.covars_[0][0]),(g2.means_[1][0],g2.covars_[1][0]))
            if aic[1] < aic[0]:
               if weight[0]*numTR > 1.0+epsilon and weight[1]*numTR > 1.0+epsilon:
                  flag = False
                  break
               else:
                  if weight[0]*numTR <= 1.0+epsilon:
                     TR_list,TR_listLabel = removeLabel(TR_list,label,0,len(TR_listNA)) #remove label 0 and return the TR list
                  elif weight[1]*numTR <= 1.0+epsilon:
                     TR_list,TR_listLabel = removeLabel(TR_list,label,1,len(TR_listNA)) #remove label 1 and return the TR list
                  TR_listNA.append(TR_listLabel)
                  #print "TR_listNA:\n",TR_listNA
            else:
               break
    tempLabelD = decouple(TR_listNA,label)
    #print "g1_stats: ",g1_stats, " g2_stats: ",g2_stats, " bic: ",bic
    #print "tempLabelD :", tempLabelD
    A2stats = andersonDarling(TR_Data[:,0])
    if aic[1] < aic[0]:
       lb = getCseparationGauss(g2_stats,stdParam)
       if lb[0] > lb[1]:
          oneAllele = False
          sM = getSimMatrix(tempLabelD)
          labelD = tempLabelD
          cl0 = [lab for lab in tempLabelD if lab == 0]
          cl1 = [lab for lab in tempLabelD if lab == 1]
          weight = [len(cl0)/float(len(cl0)+len(cl1)),len(cl1)/float(len(cl0)+len(cl1))]
          statsCall = [g2_stats[0][0],g2_stats[1][0]]
       else:
          oneAllele = True
          sM = np.zeros((1),dtype=np.int)
          weight = [1.0,0.0]
          labelD = [0 for lab in tempLabelD] 
          statsCall = [g1_stats[0],0.0]
    else:
       #print "homoz"
       lb = [0.0,0.0]
       oneAllele = True
       weight = [1.0,0.0]
       sM = np.zeros((1),dtype=np.int)
       labelD = [0 for lab in tempLabelD]
       statsCall = [g1_stats[0],0.0]
    return oneAllele,weight,labelD,sM,A2stats,lb,statsCall


def writeReadsToFile(hits,outfile,labelD,bic,numTRs,wt,stats12,numAllele,numTRna):
    index = 0
    for hit in hits:
        cl = labelD[index]
        if cl != -1:
           m = stats12[cl][0]
           stD = math.sqrt(stats12[cl][1])
           siZe = wt[cl]
           outfile.write("%s %d %d %d %d %d %.2f %s %.2f %.2f %d %d %.2f %.2f %.2f %d %d %.2f %d %d\n" %(hit.Query_id,hit.targetStart,hit.targetEnd,hit.ReftRStart,hit.ReftREnd,hit.RefPeriod,hit.RefNum,hit.tRSeq,hit.naiveTR,hit.naiveTR/hit.RefNum,numAllele,cl,siZe,m,stD,numTRs,numTRna,bic[cl],hit.anchorPre,hit.anchorSuf))
        else:
           noCluster = 'na'
           outfile.write("%s %d %d %d %d %d %.2f %s %.2f %.2f %s %d %d\n" %(hit.Query_id,hit.targetStart,hit.targetEnd,hit.ReftRStart,hit.ReftREnd,hit.RefPeriod,hit.RefNum,hit.tRSeq,hit.naiveTR,hit.naiveTR/hit.RefNum,noCluster,hit.anchorPre,hit.anchorSuf))
        index += 1
    outfile.flush()

