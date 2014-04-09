#!/usr/bin/env python
import numpy as np
import math
cimport numpy as np
cimport cython

np.import_array()

#Global Variables############
cdef int TableLen = 16000
cdef int scale = 1000
cdef double lnExpTable[16000]
############################

cdef extern from "math.h" nogil:
     long double exp(long double)
     long double log(long double) 
     int isnan(long double)
     float INFINITY

@cython.boundscheck(False)
@cython.wraparound(False)
def forwardProbLog(np.ndarray[np.int_t,ndim=1] q, int nq, np.ndarray[np.int_t,ndim=1] r, int nr, np.ndarray[np.float_t,ndim=2] Lf_M, np.ndarray[np.float_t,ndim=2] Lf_X, np.ndarray[np.float_t,ndim=2] Lf_Y, np.ndarray[np.float_t,ndim=2] eProb_M, np.ndarray[np.float_t,ndim=2] eProb_X, np.ndarray[np.float_t,ndim=1] eProb_Y, np.ndarray[np.float_t,ndim=2] tProb):

    #Variable Declarations
    cdef int i, j, qi, qi_1, rj
    cdef double emmM, emmX, emmY, LemmM, LemmX, LemmY
    cdef double LtProbMM,LtProbMX,LtProbMY,LtProbXM,LtProbXX,LtProbXY,LtProbYM,LtProbYX,LtProbYY 

    #Initialization. The Global lookup table should be initialized here...
    InitLookupTable() 

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

    #Recurrsion
    for i in range(1,nq+1):
        qi = getindicies(q[i-1]) #seq in q at i position concantenated with a gap
        if i == 1:
           emmX = 0.25 #equal prob for any base to come first. no conditional dependence for first base.
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

####Utilities for Log space####
#1.taking log
cdef inline long double eln(long double x):
    if x > 0.0:
       return log(x)
    else:
       return -1.0*INFINITY

#2.taking exp
cdef inline long double eexp(long double x):
    if isnan(x) != True:
       return exp(x)
    else:
       #print "x is nan"
       return 0.0

#3. log_sum
cdef inline long double log_sum(long double a,long double b):
    global lnExpTable
    cdef long double max, min
    cdef float scale = 1000
    max = a if (a >= b) else b
    min = b if (b < a) else a
    if (min == -1.0*INFINITY or (max-min) >= 15.7):
        return max
    else:        
        return (max + lnExpTable[int((max-min)*scale)])
    #return max if (min == -1.0*INFINITY or (max-min) >= 15.7) else (max + lnExpTable[int((max-min)*scale)])

#4. Init lookup table
cdef InitLookupTable():
    global TableLen
    global scale
    global lnExpTable
    cdef int i
    for i in range(0,TableLen):
        lnExpTable[i] = eln(1.0+eexp(float(-1.0*i/scale)))
####

cdef inline int getindicies(int q):
    cdef int qi
    if q == 65: #its an A
       qi = 0
    if q == 67: #its a C
       qi = 1
    if q == 84: #its a T
       qi = 2
    if q == 71: #its a G
       qi = 3
    return qi

@cython.boundscheck(False)
@cython.wraparound(False)
def numTRlog(np.ndarray[np.float_t,ndim=2] Lf_M, np.ndarray[np.float_t,ndim=2] Lf_X, np.ndarray[np.float_t,ndim=2] Lf_Y, int ntr, int nr, int nq):
    cdef np.ndarray[np.double_t,ndim=1] LendAt = np.zeros((nr+1))
    cdef np.ndarray[np.double_t,ndim=1] edf = np.zeros((nr+1))
    cdef int j
    cdef long double SumAllW = 0.0
    cdef long double num = 0.0
    cdef long double prob = 0.0
    cdef long double temp = -1.0*INFINITY
    cdef long double LMaxfwdProb = 0.0
    cdef long double LProbRandom = 0.0
    cdef int index = 0
    cdef long double llr = 1.0
    InitLookupTable()
    
    for j in range(1,nr+1):
        LendAt[j] = log_sum(log_sum(Lf_M[nq,j],Lf_X[nq,j]),Lf_Y[nq,j])
        if LendAt[j] >= temp:
           temp = LendAt[j]
           index = j
    #print "index float: ", <float>index
    LMaxfwdProb = temp
    LProbRandom = randomModel(nq,index)
    llr = LMaxfwdProb - LProbRandom
    if (llr > 0.0):
       SumAllW = calLogSumSeries(LendAt,nr,ntr,0)
       LogNumTR = calLogSumSeries(LendAt,nr,ntr,1) - SumAllW
       #print "LogNumTR: ",LogNumTR
       num = eexp(LogNumTR)
       #print num
       for j in range(1,nr+1):
           prob = eexp(LendAt[j]-SumAllW)
           edf[j] = prob
       #print "edf: ", edf
    else:
       num = 0.0
    return num, llr, index, LMaxfwdProb, edf, LendAt  

def getAllFwdProb(np.ndarray[np.float_t,ndim=2] Lf_M, np.ndarray[np.float_t,ndim=2] Lf_X, np.ndarray[np.float_t,ndim=2] Lf_Y, int ntr, int nr, int nq):
    InitLookupTable()
    cdef np.ndarray[np.double_t,ndim=2] LallQR = np.zeros((nq+1,nr+1))
    cdef np.ndarray[np.double_t,ndim=1] index = np.zeros((nq+1))
    cdef int i 
    cdef int j
    cdef long double temp = -1.0*INFINITY
    for i in range(1,nq+1): #The summation is performed for all i and j in nq*nr matrix:
        temp = -1.0*INFINITY
        for j in range(1,nr+1):
            LallQR[i,j] = log_sum(log_sum(Lf_M[i,j],Lf_X[i,j]),Lf_Y[i,j])
            if LallQR[i,j] >= temp: #At every i, store j where max Prob occured in index
               temp = LallQR[i,j]
               index[i] = j
    return LallQR, index

cdef long double calLogSumSeries(np.ndarray[np.double_t,ndim=1] LendAt, int nr, int ntr, int flag):
    cdef long double LogNumInit = -1.0*INFINITY
    cdef long double Lognum = 0.0
    cdef int j
    if flag == 1:
       for j in range(1,nr+1):
           Lognum = log_sum(LogNumInit,LendAt[j]+eln(j/<float>ntr))
           LogNumInit = Lognum
    else:
       for j in range(1,nr+1):
           Lognum = log_sum(LogNumInit,LendAt[j])
           LogNumInit = Lognum
    return Lognum

cdef inline long double randomModel(int nq,int nr):
    cdef long double eta = 0.05
    cdef long double g = 0.25 #prob of emitting a base = 0.25 for each base in query sequence.
    cdef long double pqr = 1.0
    cdef long double LprobRand = 1.0
    LprobRand = 2.0*log(eta) + (nq+nr+1)*log(1.0-eta) + (nq+nr)*log(g)
    return LprobRand
