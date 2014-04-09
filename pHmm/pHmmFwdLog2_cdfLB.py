#!/usr/bin/env python
import sys
import os
import numpy as np
#import pairHmmFWD2 as pHf2
import pairHmmFWD2_allNTR as pHf2
import math
import bisect
#import matplotlib.pyplot as plt


#GLOBAL VARIABLE#################################
TableLen = 16000
scale = 1000
lnExpTable = np.zeros((TableLen),dtype=np.double)
#################################################

def forwardProb(q,nq,r,nr,f_M,f_X,f_Y,eProb_M,eProb_X,eProb_Y,tProb):
    #Initialization of the probability matrices
    f_M[0,0], f_X[0,0], f_Y[0,0] = tProb[3,0], tProb[3,1], tProb[3,2]
    #f_M[0,0], f_X[0,0], f_Y[0,0] = 1.0,0.0,0.0
    for i in range(1,nq+1):
        qi = getindicies(q[i-1]) #seq in q at i position concantenated with a gap
        if i == 1:
           emm = 0.25 #a special case as emmisions in X state are conditional : RETHINK about this!!!
           f_X[i,0] = emm*(f_X[i-1,0])
        else:
           qi_1 = getindicies(q[i-2])
           emm = eProb_X[qi,qi_1]
           f_X[i,0] = emm*(tProb[0,1]*f_M[i-1,0] + tProb[1,1]*f_X[i-1,0] + tProb[2,1]*f_Y[i-1,0])
    for j in range(1,nr+1):
        rj = getindicies(r[j-1]) #seq in r at j position concantenated with a gap
        if j == 1:
           f_Y[0,j] = eProb_Y[rj]*(f_Y[0,j-1])
        else:
           f_Y[0,j] = eProb_Y[rj]*(tProb[0,2]*f_M[0,j-1] + tProb[2,2]*f_Y[0,j-1] + tProb[1,2]*f_X[0,j-1])
    for j in range(1,nr+1):
        for i in range(1,nq+1):
            qi = getindicies(q[i-1])
            rj = getindicies(r[j-1])
            if i == 1:
               emm = 0.25
            else:
               qi_1 = getindicies(q[i-2])
               emm = eProb_X[qi,qi_1]
            if i == 1 and j == 1:
               f_M[i,j] = eProb_M[qi,rj]*(f_M[i-1,j-1])
            else:
               f_M[i,j] = eProb_M[qi,rj]*(tProb[0,0]*f_M[i-1,j-1] + tProb[1,0]*f_X[i-1,j-1] + tProb[2,0]*f_Y[i-1,j-1])
            f_X[i,j] = emm*(tProb[0,1]*f_M[i-1,j] + tProb[1,1]*f_X[i-1,j] + tProb[2,1]*f_Y[i-1,j])
            f_Y[i,j] = eProb_Y[rj]*(tProb[0,2]*f_M[i,j-1] + tProb[2,2]*f_Y[i,j-1] + tProb[1,2]*f_X[i,j-1])
    return f_M,f_X,f_Y
#++++

def forwardProbLog(q,nq,r,nr,Lf_M,Lf_X,Lf_Y,eProb_M,eProb_X,eProb_Y,tProb):
    LtProbMM = eln(tProb[0,0])
    LtProbMX = eln(tProb[0,1])
    LtProbMY = eln(tProb[0,2])

    LtProbXM = eln(tProb[1,0])
    LtProbXX = eln(tProb[1,1])
    LtProbXY = eln(tProb[1,2])
    
    LtProbYM = eln(tProb[2,0])
    LtProbYX = eln(tProb[2,1])
    LtProbYY = eln(tProb[2,2])
   
    Lf_M[0,0], Lf_X[0,0], Lf_Y[0,0] = eln(tProb[3,0]), eln(tProb[3,1]), eln(tProb[3,2])

    for i in range(1,nq+1):
        qi = getindicies(q[i-1]) #seq in q at i position concantenated with a gap
        if i == 1:
           emmX = 0.25
           Lf_X[i,0] = eln(emmX) + Lf_X[i-1,0]
        else:
           qi_1 = getindicies(q[i-2])
           emmX = eProb_X[qi,qi_1]
           LemmX = eln(emmX)
           Lf_X[i,0] = LemmX + log_sum(log_sum(LtProbMX+Lf_M[i-1,0],LtProbXX+Lf_X[i-1,0]),LtProbYX+Lf_Y[i-1,0]) 

    for j in range(1,nr+1):
        rj = getindicies(r[j-1]) #seq in r at j position concantenated with a gap
        LemmY = eln(eProb_Y[rj])
        if j == 1:
           Lf_Y[0,j] = LemmY + Lf_Y[0,j-1]
        else:
           Lf_Y[0,j] = LemmY + log_sum(log_sum(LtProbMY+Lf_M[0,j-1],LtProbYY+Lf_Y[0,j-1]),LtProbXY+Lf_X[0,j-1]) 

    for j in range(1,nr+1):
        for i in range(1,nq+1):
            qi = getindicies(q[i-1])
            rj = getindicies(r[j-1])
            LemmM = eln(eProb_M[qi,rj])
            LemmY = eln(eProb_Y[rj])
            if i == 1:
               LemmX = eln(0.25)
            else:
               qi_1 = getindicies(q[i-2])
               LemmX = eln(eProb_X[qi,qi_1])
            if i == 1 and j == 1:
               Lf_M[i,j] = LemmM + Lf_M[i-1,j-1]
            else:
               Lf_M[i,j] = LemmM + log_sum(log_sum(LtProbMM+Lf_M[i-1,j-1],LtProbXM+Lf_X[i-1,j-1]),LtProbYM+Lf_Y[i-1,j-1])
            Lf_X[i,j] = LemmX + log_sum(log_sum(LtProbMX+Lf_M[i-1,j],LtProbXX+Lf_X[i-1,j]),LtProbYX+Lf_Y[i-1,j])
            Lf_Y[i,j] = LemmY + log_sum(log_sum(LtProbMY+Lf_M[i,j-1],LtProbXY+Lf_X[i,j-1]),LtProbYY+Lf_Y[i,j-1])
    return Lf_M,Lf_X,Lf_Y

###UTILITY FUNCTIONS FOR LOG/EXP SPACE###

def log_sum_prev(a,b):
    if a > b:
       return a + eln(1.0 + eexp(b-a))
    else:
       return b + eln(1.0 + eexp(a-b))

def InitLookupTable():
    for i in range(0,TableLen):
        lnExpTable[i] = eln(1.0+eexp((float(-1.0*i/scale))))
    #print "Message: LookUp Table initiated"
#def LogExpLookup(x):
    

def log_sum(a,b):
    max = a if (a >= b) else b
    min = b if (b < a) else a
    #return max if (min == -1.0*float("inf") or (max-min) >= 15.7) else (max + eln(1.0 + eexp(min-max)))
    return max if (min == -1.0*float("inf") or (max-min) >= 15.7) else (max + lnExpTable[int((max-min)*scale)])

def eexp(x):
    if math.isnan(x) != True:
       return np.exp(x)
    else:
       print "x is nan"
       return 0

def eln(x):
    if x > 0.0:
       return np.log(x)
    else:
       print "x <= 0.0 encountered"
       return -1.0*float("inf")

########################################
def getindicies(q):
    if q == 65: #its an A
       qi = 0
    if q == 67: #its a C
       qi = 1
    if q == 84: #its a T
       qi = 2
    if q == 71: #its a G
       qi = 3
    return qi


def numTRlog(Lf_M,Lf_X,Lf_Y,ntr,nr,nq):
    LendAt = np.zeros((nr+1),dtype=np.float)
    endAt = np.zeros((nr+1),dtype=np.float)
    temp = -1.0*float("inf")
    index = 0
    edf = []
    sumEndAt = 0.0
    for j in range(1,nr+1):
        #if j%ntr == 0:
           endAt[j] = (np.exp(Lf_M[nq,j])+np.exp(Lf_X[nq,j])+np.exp(Lf_Y[nq,j])) #this gives right answer...
           sumEndAt += endAt[j]
           LendAt[j] = log_sum(log_sum(Lf_M[nq,j],Lf_X[nq,j]),Lf_Y[nq,j])
           if LendAt[j] >= temp:
              temp = LendAt[j]
              index = j
    print "LendAt: ", LendAt
    print "endAt: ", endAt
    LMaxfwdProb = temp
    LProbRandom = randomModel(nq,index)
    print "Max forward probability in Log space: ", LMaxfwdProb
    print "Occurs at index: ",index
    llr = LMaxfwdProb - LProbRandom
    print "Log-Odds ratio Score: ", llr

    if (llr > 0.0):
       expTR = 0.0
       for j in range(1,nr+1):
           expTR += j*endAt[j]
       expTR = expTR/sumEndAt
       print "expTR: ", (expTR/float(ntr)) 
       SumAllW = calLogSumSeries(LendAt,nr,ntr,0)
       print "All W sum in Log: ",(SumAllW)#, "sumL: ",sumL
       LogNumTR = calLogSumSeries(LendAt,nr,ntr,1) - SumAllW
       num = eexp(LogNumTR)
       print "num: ", num
       for j in range(1,nr+1):
           #if j%ntr == 0:
              prob = eexp(LendAt[j]-SumAllW)
              #if prob > 0.05: #significance level of 5% in a way...(think about this threshold)!
              edf.append(prob)
    else:
       num = 0.0
       print "Log-Odds ratio Score < 0.0"
    return num, llr, index, LProbRandom, edf,LendAt

def generateSamples(edfArray,index,ntr,num):
    samples, edf = [], []
    edf1 = edfArray.tolist()
    for i in range(0,len(edf1)):
        if edf1[i] > 0.01:
           edf.append((edf1[i],i))
    #print "edf...: ", edf
    if len(edf) != 0:
       cdf_Edf, cdf = [], []
       sumF = 0.0
       intervalInit = edf[0][1]-1
       for wj_k in edf:
           sumF += wj_k[0]
           cdf_Edf.append((sumF,(intervalInit,wj_k[1])))
           cdf.append(sumF)
           intervalInit = wj_k[1]
       #print "cdf_Edf: ",cdf_Edf
       for i in range(num):
           temp = np.random.random_sample() #results are from the continuous uniform distribution
           #print "temp: ",temp
           ind = bisect.bisect_left(cdf_Edf,temp)
           ind = bisect.bisect_left(cdf,temp)
           #print "ind: ",ind
           if ind == len(cdf):
              interval = (cdf_Edf[ind-1][1][1],cdf_Edf[ind-1][1][1]+1)
           else:
              interval = cdf_Edf[ind][1]
           #print "interval: ",interval
           samples.append(((interval[0]+interval[1])/2.0)/float(ntr))
       #print "samples: ",samples
    return samples

def getAllFwdProb(Lf_M,Lf_X,Lf_Y,ntr,nr,nq):
    LallQR = np.zeros((nq+1,nr+1),dtype=np.float)
    index = np.zeros((nq+1),dtype=np.float)
    for i in range(1,nq+1): #The summation is performed for all i and j in nq*nr matrix:
        temp = -1.0*float("inf")
        for j in range(1,nr+1):
            if j%ntr == 0:
               LallQR[i,j] = log_sum(log_sum(Lf_M[i,j],Lf_X[i,j]),Lf_Y[i,j]) 
               if LallQR[i,j] >= temp: #At every i, store j where max Prob occured in index
                  temp = LallQR[i,j]
                  index[i] = j
    return LallQR,index

def unionRegion12(region1,region2):
    if len(region1) != 0:
       merged12 = []
       mergedHighDiv = []
       merged12 = list(set(region1).union(set(region2)))
       print "merged12: ",merged12
       merged12.sort()
       print "merged12: ",merged12
       idx_list = []
       for i in range(0,len(region2)):
           idx = merged12.index(region2[i])
           if idx == 0: 
              print "idx: ",idx
              if merged12[idx+1][0] == region2[i][1] + 1:
                 idx_list.append(i)
                 mergedHighDiv.append((region2[i][0],merged12[idx+1][1]))
           else:
              print "idx: ",idx
              if idx == len(merged12):
                  if merged12[idx-1][1] == region2[i][0] - 1:
                     idx_list.append(i)
                     mergedHighDiv.append((merged12[idx-1][0],region2[i][1]))
              else:
                  if merged12[idx+1][0] == region2[i][1] + 1:
                     idx_list.append(i)
                     if merged12[idx-1][1] == region2[i][0] - 1:
                        mergedHighDiv.append((merged12[idx-1][0],merged12[idx+1][1]))
                     else:
                        mergedHighDiv.append((region2[i][0],merged12[idx+1][1]))
                  if merged12[idx-1][1] == region2[i][0] - 1:
                     idx_list.append(i)
                     mergedHighDiv.append((merged12[idx-1][0],region2[i][1]))
           print "idx_list: ",idx_list
       print "merged High: ", mergedHighDiv
       print "idx_list: ",idx_list
       mergedHighFinal = []
       for i in range(0,len(idx_list)):
           idx = idx_list[i]
           region2[idx] = mergedHighDiv[i]
       print "region2: ",region2
    return region2

def getLowerBound2(LallQR,index,nq,ntr,LLRcutOff):
    region1,region2 = getDiffLLR2(LallQR,index,nq,ntr,LLRcutOff)
    #print "region2: ",region2
    lowerBound = 0.0
    highDivRegion = []
    #mergedHighDiv = []
    if len(region2) != 0:
       highDivRegion = consolidateRegion(region2)
       #lowDivRegion = consolidateRegion(region1)
       #print "region1: ", lowDivRegion
       #print "highDivRegion: ", highDivRegion
       #mergedHighDiv = unionRegion12(lowDivRegion,highDivRegion)
       #print "mergedHighDiv: ",mergedHighDiv
       #mergedHighDiv = consolidateRegion(mergedHighDiv)
       #print "mergedHighDiv: ",mergedHighDiv
       #highDivLen = regionLen(mergedHighDiv)
       highDivLen = regionLen(highDivRegion)
       #print "highDivLen: ",highDivLen
       lowerBound = highDivLen/float(ntr)
       #print "lowerBound (ntr): ",lowerBound
    return lowerBound, highDivRegion

#def getLowerBound2(LallQR,index,nq,ntr,LLRcutOff):
#    region1,region2 = getDiffLLR2(LallQR,index,nq,ntr,LLRcutOff)
#    lowerBound = 0.0
#    print "Have reg1 reg2"
#    highDivRegion = []
#    if len(region2) != 0:
#       highDivRegion = consolidateRegion(region2)
#       #print "highDivRegion: ", highDivRegion
#       highDivLen = regionLen(highDivRegion)
#       lowerBound = highDivLen/ntr
#       print "length of non-TR sequences mod ntr: ", lowerBound
#    return lowerBound, highDivRegion

def regionLen(region):
    highDivLen = 0.0
    for i in range(0,len(region)):
        highDivLen += region[i][1] - region[i][0] + 1
    return highDivLen

def consolidateRegion(region):
    contRegion = []
    if len(region) != 0:
       randomRegion = region[0]
       for i in range(0,len(region)):
           if i != 0:
              if randomRegion[1] + 1 == region[i][0] or randomRegion[1] == region[i][0]:
                 randomRegion = (randomRegion[0],region[i][1])
              else:
                 contRegion.append(randomRegion)
                 randomRegion = region[i]
       contRegion.append(randomRegion)
    return contRegion


def getDiffLLR_OldMethod(LallQR,index,nq,ntr):
    temp = -1.0*float("inf")
    llr = np.zeros((nq+1),dtype=np.float)
    diffLLR = np.zeros((nq+1),dtype=np.float)
    plt_d, plt_x, plt_llr = [], [], []
    for i in range(1,nq+1):
        idx = index[i]  #get j where max forward prob occurs for a given i
        randMod = randomModel(i,idx)
        llr[i] = LallQR[i,idx] - randMod  #get the difference between Max fwd Prob and prob from random model here
        diffLLR[i] = llr[i] - llr[i-1]
        plt_x.append(i)
        plt_d.append(diffLLR[i])
        plt_llr.append(llr[i])
        if llr[i] > 0:
           if llr[i] >= temp:
                     temp = llr[i]
                     maxLlr = i
    return diffLLR, llr, plt_d, plt_x, plt_llr

def getDiffLLR2(LallQR,index,nq,ntr,LLRcutOff):
    temp = -1.0*float("inf")
    llr = np.zeros((nq+1),dtype=np.float)
    idx_prev = index[0]
    idx_prev_i = 0
    region1,region2 = [],[]
    for i in range(1,nq+1):
        idx = index[i]  #get j where max forward prob occurs for a given i
        randMod = randomModel(i,idx)
        llr[i] = LallQR[i,idx] - randMod  #get the difference between Max fwd Prob and prob from random model here
        if idx > idx_prev:
           diffLLR_idx = llr[i] - llr[idx_prev_i]
           diff_i = i - idx_prev_i
           #print "diff w.r.t idx: ",diffLLR_idx," for idx ",idx," and idx_prev ",idx_prev," at i: ",i, " llrIndex at i: ",llr[i]," diff in i: ",diff_i," derivativeLLR: ",(diffLLR_idx/diff_i)
           if diff_i > ntr:
              if diffLLR_idx >= LLRcutOff:
                 region1.append((idx_prev_i+1,i)) #append Region1
              else:
                 region2.append((idx_prev_i+1,i)) #append Region2
           idx_prev = idx
           idx_prev_i = i
    lowDivRegion = consolidateRegion(region1)
    #if len(lowDivRegion) != 0:
    #   #print "Low divergence region (somewhat error prone) Length: ", regionLen(region1)
    return region1,region2

def makePlotsLLR(x,d,ll,qName):
    plt.plot(x,d)
    plt.grid(True)
    l = plt.axhline(linewidth=1,color='r')
    #plt.plot(x,ll)
    #   plt.axvline(x=Idx,linewidth=2,color='k')
    #outDir = "/projects/HuPac/repartition_MAR29/simulation/SVrun3/newSv"
    #fN = qName + "_" + "querySignal.png"
    #fO = outDir + "/" + fN
    #plt.savefig(fO)
    plt.show()

def getLowerBound_OldMethod(LallQR,index,nq,ntr,LLRcutOff,qName):
    diffLLR,llr,plt_d,plt_x,plt_llr = getDiffLLR_OldMethod(LallQR,index,nq,ntr)
    #makePlotsLLR(plt_x,plt_d,plt_llr,qName)
    return diffLLR,llr 

def getLowerBound3(LallQR,index,nq,ntr,LLRcutOff,qName):
    region1,region2,diffLLR,llr,plt_d,plt_x,plt_llr = getDiffLLR3(LallQR,index,nq,ntr,LLRcutOff)
    #makePlotsLLR(plt_x,plt_d,plt_llr,qName)
    #print "region2: ",region2
    lowerBound = 0.0
    highDivRegion = []
    #mergedHighDiv = []
    if len(region2) != 0:
       highDivRegion = consolidateRegion(region2)
       #lowDivRegion = consolidateRegion(region1)
       #print "region1: ", lowDivRegion
       #print "highDivRegion: ", highDivRegion
       #mergedHighDiv = unionRegion12(lowDivRegion,highDivRegion)
       #print "mergedHighDiv: ",mergedHighDiv
       #mergedHighDiv = consolidateRegion(mergedHighDiv)
       #print "mergedHighDiv: ",mergedHighDiv
       #highDivLen = regionLen(mergedHighDiv)
       highDivLen = regionLen(highDivRegion)
       #print "highDivLen: ",highDivLen
       lowerBound = highDivLen/float(ntr)
       #print "lowerBound (ntr): ",lowerBound
    return lowerBound, highDivRegion,diffLLR,llr

def getDiffLLR3(LallQR,index,nq,ntr,LLRcutOff):
    temp = -1.0*float("inf")
    llr = np.zeros((nq+1),dtype=np.float)
    diffLLR = np.zeros((nq+1),dtype=np.float)
    idx_prev = index[0]
    idx_prev_i = 0
    region1,region2 = [],[]
    plt_d, plt_x, plt_llr = [], [], []
    for i in range(1,nq+1):
        idx = index[i]  #get j where max forward prob occurs for a given i
        randMod = randomModel(i,idx)
        llr[i] = LallQR[i,idx] - randMod  #get the difference between Max fwd Prob and prob from random model here
        if idx%ntr == 0: #check to see if the j where max P(O|lambda,j) occurs is a multiple of ntr (the length of the TR element), if not, then that means that j* has not incremented by ntr
           diffLLR_idx = llr[i] - llr[idx_prev_i] #if j increases, then we need to calculate the diff in LLR from previous value and this value
           diffLLR[i] = diffLLR_idx
           diff_i = i - idx_prev_i
           plt_x.append(i)
           plt_d.append(diffLLR[i])
           plt_llr.append(llr[i])
           if diffLLR_idx < LLRcutOff: #if the rate of change of LLR is less than cutoff, that implies a non-TR region, with possible SVs
              region2.append((idx_prev_i,i)) #append Region2
              #print "region2", region2
           else:
              if diff_i > ntr: #Even if the rate of change of LLR is more than cutoff, but this change occurs when i increments by more than ntr, it implies error prone-region or region1. This coule be a putative SV region.
                 region1.append((idx_prev_i,i)) #append Region1
                 #print "region1", region1
           idx_prev = idx
           idx_prev_i = i
    return region1,region2,diffLLR,llr,plt_d,plt_x,plt_llr

def getDiffLLR3b(LallQR,index,nq,ntr,LLRcutOff):
    temp = -1.0*float("inf")
    llr = np.zeros((nq+1),dtype=np.float)
    diffLLR = np.zeros((nq+1),dtype=np.float)
    print "length of index: ", index.size
    print index
    #generate index mod ntr
    idx_prev = index[0]
    idx_prev_i = 0
    region1,region2 = [],[]
    plt_d, plt_x, plt_llr = [], [], []
    for i in range(1,nq+1):
        idx = index[i]  #get j where max forward prob occurs for a given i
        randMod = randomModel(i,idx)
        llr[i] = LallQR[i,idx] - randMod  #get the difference between Max fwd Prob and prob from random model here
        print "i: ", i, " idx: ", idx
        if idx > idx_prev: #check if the 'j' or number of TRs index increase from previous
           if idx%ntr == 0:
              diffLLR_idx = llr[i] - llr[idx_prev_i] #if j increases, then we need to calculate the diff in LLR from previous value and this value
              diffLLR[i] = diffLLR_idx
              diff_i = i - idx_prev_i
              print "diff w.r.t idx: ",diffLLR_idx," diff in i: ",diff_i," i: ",i," idx_prev_i: ",idx_prev_i
              plt_x.append(i)
              plt_d.append(diffLLR[i])
              plt_llr.append(llr[i])
              if diff_i > ntr:
                 print "I am here: diff_i > ntr"
                 print "diff w.r.t idx: ",diffLLR_idx," diff in i: ",diff_i," i: ",i," idx_prev_i: ",idx_prev_i
                 if diffLLR_idx >= LLRcutOff:
                    region1.append((idx_prev_i+1,i)) #append Region1
                    print "region1", region1
                 else:
                    region2.append((idx_prev_i+1,i)) #append Region2
                    print "region2"
              idx_prev = idx
              idx_prev_i = i
        if llr[i] > 0:
           if llr[i] >= temp:
                     temp = llr[i]
                     maxLlr = i
    lowDivRegion = consolidateRegion(region1)
    #if len(lowDivRegion) != 0:
    #   #print "Low divergence region (somewhat error prone) Length: ", regionLen(region1)
    return region1,region2,diffLLR,llr,plt_d,plt_x,plt_llr,llr[maxLlr],maxLlr

def calLogSumSeries(LendAt,nr,ntr,flag):
    LogNumInit = -1.0*float("inf")
    Lognum = 0.0
    if flag == 1:
       for j in range(1,nr+1):
           #if j%ntr == 0:
              Lognum = log_sum(LogNumInit,LendAt[j]+eln(j/float(ntr)))
              LogNumInit = Lognum
    else:
       for j in range(1,nr+1):
           #if j%ntr == 0:
              Lognum = log_sum(LogNumInit,LendAt[j])
              LogNumInit = Lognum
    return Lognum

def randomModel(nq,nr):
    eta = 0.05 #Some small parameter ensures that first sequence (q in our case) is emitted by the first state. But, also need to understand how to select this variable more tightly.
    g = 0.25 #prob of emitting a base = 0.25 (from alphabet set = {A,C,G,T}) for each base in query sequence.
    pqr = 1.0
    LprobRand = 2.0*np.log(eta) + (nq+nr+1)*np.log(1.0-eta) + (nq+nr)*np.log(g)
    #LprobRand = 2.0*np.log(eta) + (nq+nr+1)*np.log(1.0-eta) + (nq)*np.log(g)
    #print "LprobRand: ",LprobRand
    return LprobRand

#@profile
def pairHmmRun(QuerySeq,trSeq,tProb,eProb_M,eProb_X):
    ###
    InitLookupTable()
 
    #print "QuerySeq ", QuerySeq
    #print "trSeq ", trSeq
    #print "tProb:\n", tProb
    pen = 1.0
    eProb_Y = np.array([pen, pen, pen, pen])

    nq = len(QuerySeq) #QuerySeq represents the putative TR region in the Query from 3DP step
    ntr = len(trSeq) #trSeq is the repeat Element

    #coarse estimation of TRs in QuerySeq
    #estTRinQ = int(nq/ntr)
    estTRinQ = (nq/float(ntr))
    #A factor reflecting unrestricted TR estimation
    mltFactor = 1.5

    #Converting QuerySeq to integers
    q = np.arange(nq,dtype=np.int)
    tempQ = map(ord,QuerySeq)
    q = np.asarray(tempQ)

    #Defining the repeatSeq for pairHmm
    repeatSeq = trSeq*(int(mltFactor*estTRinQ)+1)
    #print "estTRinQ: ",estTRinQ*mltFactor
    #print "repeatSeq: ",repeatSeq 
    #Constructing the repeat array and converting to integers
    nr = len(repeatSeq)
    r = np.arange(nr,dtype=np.int)
    tempR = map(ord,repeatSeq)
    r = np.asarray(tempR)

    #print "nr: ",nr," nq: ",nq, " ntr: ", ntr   

    #np.set_printoptions(suppress=True)
    Lf_M = np.zeros((nq+1,nr+1))-1.0*float("inf")
    Lf_X= np.zeros((nq+1,nr+1))-1.0*float("inf")
    Lf_Y = np.zeros((nq+1,nr+1))-1.0*float("inf")
    Lf_M,Lf_X,Lf_Y = pHf2.forwardProbLog(q,nq,r,nr,Lf_M,Lf_X,Lf_Y,eProb_M,eProb_X,eProb_Y,tProb)
    print "##done_FWD_Calc##\n"

    edf = np.zeros((nr+1))
    num_TR,llr,index,ProbMax,edf,LendAt = pHf2.numTRlog(Lf_M,Lf_X,Lf_Y,ntr,nr,nq)
    #num_TR,llr,index,ProbMax,edf,LendAt = numTRlog(Lf_M,Lf_X,Lf_Y,ntr,nr,nq)
    #print "LendAt: ", LendAt
    #print "index: ", index
    #print "edf: ", edf
    samplesNumTR = generateSamples(edf,index,ntr,50) 
    print "##sampling##\n"

    LLRcutOff = 0.01
    LallQR,index = pHf2.getAllFwdProb(Lf_M,Lf_X,Lf_Y,ntr,nr,nq) 
    #print "index where max prob exists: ", index
    print "##AllFwdProb##\n"
    #lowerBound, devList, maxLlrIdx = getLowerBound(LallQR,index,nq,ntr,LLRcutOff)
    highDivLen, highDivRegion = getLowerBound2(LallQR,index,nq,ntr,LLRcutOff)
    print "##LowerBound##\n"
    
    return num_TR, llr, ProbMax, samplesNumTR, highDivLen, highDivRegion

def pairHmmRun2(QuerySeq,trSeq,tProb,eProb_M,eProb_X,qName):
    ###
    InitLookupTable()

    #print "QuerySeq ", QuerySeq
    #print "trSeq ", trSeq
    #print "tProb:\n", tProb
    pen = 1.0
    eProb_Y = np.array([pen, pen, pen, pen])

    nq = len(QuerySeq) #QuerySeq represents the putative TR region in the Query from 3DP step
    ntr = len(trSeq) #trSeq is the repeat Element

    #coarse estimation of TRs in QuerySeq
    #estTRinQ = int(nq/ntr)
    estTRinQ = (nq/float(ntr))
    #A factor reflecting unrestricted TR estimation
    mltFactor = 1.5

    #Converting QuerySeq to integers
    q = np.arange(nq,dtype=np.int)
    tempQ = map(ord,QuerySeq)
    q = np.asarray(tempQ)

    #Defining the repeatSeq for pairHmm
    repeatSeq = trSeq*(int(mltFactor*estTRinQ)+1)
    #print "estTRinQ: ",estTRinQ*mltFactor
    #print "repeatSeq: ",repeatSeq 
    #Constructing the repeat array and converting to integers
    nr = len(repeatSeq)
    r = np.arange(nr,dtype=np.int)
    tempR = map(ord,repeatSeq)
    r = np.asarray(tempR)

    #print "nr: ",nr," nq: ",nq, " ntr: ", ntr   

    #np.set_printoptions(suppress=True)
    Lf_M = np.zeros((nq+1,nr+1))-1.0*float("inf")
    Lf_X= np.zeros((nq+1,nr+1))-1.0*float("inf")
    Lf_Y = np.zeros((nq+1,nr+1))-1.0*float("inf")
    Lf_M,Lf_X,Lf_Y = pHf2.forwardProbLog(q,nq,r,nr,Lf_M,Lf_X,Lf_Y,eProb_M,eProb_X,eProb_Y,tProb)
    print "##done_FWD_Calc##\n"

    edf = np.zeros((nr+1))
    num_TR,llr,index,ProbMax,edf,LendAt = pHf2.numTRlog(Lf_M,Lf_X,Lf_Y,ntr,nr,nq)

    LLRcutOff = 0.01
    LallQR,index = pHf2.getAllFwdProb(Lf_M,Lf_X,Lf_Y,ntr,nr,nq)
    print "##AllFwdProb##\n"

    #lowerBound, devList, diffLLR,LLR = getLowerBound3(LallQR,index,nq,ntr,LLRcutOff,qName)
    diffLLR,LLR = getLowerBound_OldMethod(LallQR,index,nq,ntr,LLRcutOff,qName)
    print "##LowerBound##\n"

    #return num_TR, diffLLR, LLR,lowerBound,devList #UNCOMMENT THIS FOR THE NEW METHOD
    return num_TR, diffLLR, LLR #THIS IS FOR THE OLD METHOD

def getGlobalProb(target):
   #all2_B_chrY.m5.eProbM.txt
    path2params = "/hpc/users/ummata01/gitrepos/workIn/debug/pacBioTR/psCount/parameters/"
    prefix = "all2_B_"
    suffixEM = ".m5.eProbM.txt"
    suffixEX = ".m5.eProbX.txt"
    suffixT = ".m5.tProb.txt"
    suffixI = ".m5.InitProb.txt"
    fn = path2params + prefix + target
    tProb = np.loadtxt(fn + suffixT)
    eProbM = np.loadtxt(fn + suffixEM)
    eProbX = np.loadtxt(fn + suffixEX)
    ProbInit = np.loadtxt(fn + suffixI)

    temp = np.vstack((tProb,ProbInit))
    colOfZero = np.zeros((4,1))
    temp2 = np.hstack((temp,colOfZero))
    tau = np.array([[0.0],[0.0],[0.0],[0.0]])
    temp3 = np.hstack((temp2,tau))
    rowOfZero = np.zeros((1,5))
    tProbHmm = np.vstack((temp3,rowOfZero))
    return tProbHmm, eProbM, eProbX

if  __name__== '__main__':

    inm5fn = sys.argv[1] #input SV file 

    tempm5 = inm5fn.split("/")[-1].split(".")[0]
    outBinfn = sys.argv[2] #directory to store the .binned_anchors files
    hgT_outFn = outBinfn + "/" + tempm5 + ".toHMM.sv.txt"
    outfile = open(hgT_outFn, 'w')

    tProb,eProb_M,eProb_X = getGlobalProb('chr1')

    for hline in open(inm5fn, 'r'):
        hvalues = hline.rstrip("\n").split()
        qName = hvalues[0]
        print qName
        q = hvalues[1]
        lenQ = len(q)
        tr = hvalues[2]
        num_TR, diffLLR, LLR = pairHmmRun2(q,tr,tProb,eProb_M,eProb_X,qName)
        diffLLRstr = ' '.join(map(str,diffLLR))
        LLRstr = ' '.join(map(str,LLR))
        #outfile.write("%s %d %.2f %s %s %.2f %s\n"%(qName,lenQ,num_TR,diffLLRstr,LLRstr,lowerBound,highDivRegion))
        outfile.write("%s %d %.2f %s %s\n"%(qName,lenQ,num_TR,diffLLRstr,LLRstr))
        outfile.flush()
    outfile.close()
