#!/usr/bin/env python
import os
import sys
import numpy as np
import countSeq2 as cS

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

def balancetProb(tProb): 
    for i in [0,1,2]:
        for j in [0,1,2]:
            if tProb[i,j] < 0.001:
               tProb[i,j] = 0.001
        sum_i = np.sum(tProb[i,:])
        tProb[i,0],tProb[i,1],tProb[i,2] = tProb[i,0]/sum_i,tProb[i,1]/sum_i,tProb[i,2]/sum_i
    return tProb


def balanceProbInit(tProb): 
    for i in [0,1,2]:
        if tProb[i] < 0.001:
           tProb[i] = 0.001
    sum_i = np.sum(tProb)
    tProb[0],tProb[1],tProb[2] = tProb[0]/sum_i,tProb[1]/sum_i,tProb[2]/sum_i
    return tProb

def countsOnList(SqPq):
    countsTotal = np.zeros((5,5))
    countsTotalEmm = np.zeros((4,4))
    transTotal = np.zeros((3,3))
    gammaTotal = np.zeros((3))
    for hit in SqPq:
        if hit != '':
           alnQ,alnT = hit[1],hit[0]
           nq = len(alnQ)
           q = np.arange(nq,dtype=np.int_)
           tempQ = map(ord,alnQ)
           q = np.asarray(tempQ)
           nt = len(alnT)
           t = np.arange(nt,dtype=np.int_)
           tempT = map(ord,alnT)
           t = np.asarray(tempT)
           countsLocal,transLocal,countsLocalEmm,gammaLocal = cS.getCountsFromAln(q,nq,t,nt)
           countsTotal = np.add(countsLocal,countsTotal)
           countsTotalEmm = np.add(countsLocalEmm,countsTotalEmm)
           transTotal = np.add(transLocal,transTotal)
           gammaTotal = np.add(gammaLocal,gammaTotal)
#    print "countsTotal:\n",countsTotal
           #print "countsTotalEmm:\n",countsTotalEmm

    tProb = cS.getTransProb(transTotal)
#    print "tProb:\n",tProb
    tProb = balancetProb(tProb) 
#    print "tProbB:\n",tProb
    ProbInit = cS.getInitProb(gammaTotal)
#    print "ProbInit:\n",ProbInit
    ProbInit = balanceProbInit(ProbInit) 
#    print "ProbInitB:\n",ProbInit
    temp = np.vstack((tProb,ProbInit))
    colOfZero = np.zeros((4,1))
    temp2 = np.hstack((temp,colOfZero))
    tau = np.array([[0.0],[0.0],[0.0],[0.0]])
    temp3 = np.hstack((temp2,tau))
    rowOfZero = np.zeros((1,5))
    tProbHmm = np.vstack((temp3,rowOfZero))
#    print "Transition Prob.:\n",tProbHmm

    eProbM = cS.new_eProb(countsTotal)
#    print "eProb in match state:\n",eProbM

    eProbX = cS.new_eProb(countsTotalEmm)
#    print "eProb in ins state:\n",eProbX
    return tProbHmm, eProbM, eProbX 

def countsOnAlignment(qXa):
    countsTotal = np.zeros((5,5))
    countsTotalEmm = np.zeros((4,4))
    transTotal = np.zeros((3,3))
    gammaTotal = np.zeros((3))
    alnQ,alnT = qXa[1],qXa[0] 
    nq = len(alnQ)
    q = np.arange(nq,dtype=np.int_)
    tempQ = map(ord,alnQ)
    q = np.asarray(tempQ)
    nt = len(alnT)
    t = np.arange(nt,dtype=np.int_)
    tempT = map(ord,alnT)
    t = np.asarray(tempT)
    countsLocal,transLocal,countsLocalEmm,gammaLocal = cS.getCountsFromAln(q,nq,t,nt)
    countsTotal = np.add(countsLocal,countsTotal)
    countsTotalEmm = np.add(countsLocalEmm,countsTotalEmm)
    transTotal = np.add(transLocal,transTotal)
    gammaTotal = np.add(gammaLocal,gammaTotal)

    print "countsTotal:\n",countsTotal
    tProb = cS.getTransProb(transTotal)
    print "tProb:\n",tProb
    ProbInit = cS.getInitProb(gammaTotal)
    print "InitProb:\n", ProbInit
    eProbM = cS.new_eProb(countsTotal)
    print "eProb in match state:\n",eProbM
    eProbX = cS.new_eProb(countsTotalEmm)
    print "eProb in ins state:\n",eProbX
