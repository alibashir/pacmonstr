#!/usr/bin/env python
import numpy as np
cimport numpy as np
cimport cython

np.import_array()
@cython.boundscheck(False)
@cython.wraparound(False)
def getCountsFromAln(np.ndarray[np.int_t,ndim=1,negative_indices=False] q, int nq, np.ndarray[np.int_t,ndim=1,negative_indices=False] t, int nt):
    #0=m,1=x,2=y (on both axis)
    cdef np.ndarray[np.float_t,ndim=2,negative_indices=False] trans = np.zeros((3,3),dtype=np.float)
    #0=A,1=C,2=T,3=G,4='-' (on both axis)
    cdef np.ndarray[np.float_t,ndim=2,negative_indices=False] countsLocal = np.zeros((5,5),dtype=np.float)
    cdef np.ndarray[np.float_t,ndim=2,negative_indices=False] countsLocalEmm = np.zeros((4,4),dtype=np.float)
    cdef np.ndarray[np.float_t,ndim=1,negative_indices=False] gamma = np.zeros((3),dtype=np.float)
    cdef int i = 0
    cdef int q_i = 0
    cdef int t_i = 0
    cdef int q_i1 = 0
    cdef int t_i1 = 0
    if nt == nq:
       for i in range(0,nq-1):
           q_i = getindicies(q[i])
           t_i = getindicies(t[i])
           q_i1 = getindicies(q[i+1])
           t_i1 = getindicies(t[i+1])
           if i == 0:
              if  (q_i != 4) and (t_i != 4):
                  gamma[0] += 1
              if  (q_i == 4) and (t_i != 4):
                  gamma[2] += 1
              if  (q_i != 4) and (t_i == 4):
                  gamma[1] += 1
           if  (q_i != 4) and (t_i != 4):
               if  (q_i1 != 4) and (t_i1 != 4):
                   trans[0,0] += 1
               if  (q_i1 != 4) and (t_i1 == 4):
                   trans[0,1] += 1
                   countsLocalEmm[q_i,q_i1] += 1
               if  (q_i1 == 4) and (t_i1 != 4):
                   trans[0,2] += 1
           if  (q_i == 4) and (t_i != 4):
               if  (q_i1 != 4) and (t_i1 != 4):
                   trans[2,0] += 1
               if  (q_i1 != 4) and (t_i1 == 4):
                   trans[2,1] += 1
               if  (q_i1 == 4) and (t_i1 != 4):
                   trans[2,2] += 1
           if  (q_i != 4) and (t_i == 4):
               if  (q_i1 != 4) and (t_i1 != 4):
                   trans[1,0] += 1
                   if q_i1 == q_i and q_i1 == t_i1:
                      countsLocalEmm[q_i,q_i1] += 1
               if  (q_i1 != 4) and (t_i1 == 4):
                   trans[1,1] += 1
                   countsLocalEmm[q_i,q_i1] += 1
               if  (q_i1 == 4) and (t_i1 != 4):
                   trans[1,2] += 1
           countsLocal[q_i,t_i] += 1
    return countsLocal, trans, countsLocalEmm, gamma

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
    if q == 45: #its a '-'
       qi = 4
    return qi

@cython.boundscheck(False)
@cython.wraparound(False)
def getTransProb(np.ndarray[np.float_t,ndim=2,negative_indices=False] T):
    cdef np.ndarray[np.float_t,ndim=2,negative_indices=False] tProb = np.zeros((3,3),dtype=np.float)
    cdef np.ndarray[np.float_t,ndim=1,negative_indices=False] sumT = np.zeros((3),dtype=np.float)
    cdef int i = 0
    cdef int j = 0
    cdef float epsilon = 0.001
    for i in range(0,3):
        sumT[i] = T[i,0]+T[i,1]+T[i,2]
        if sumT[i] != 0.0:
           for j in range(0,3):
               tProb[i,j] = T[i,j]/(sumT[i])
        else:
           tProb[i,0], tProb[i,1], tProb[i,2] = 1.0-2.0*epsilon,epsilon,epsilon
    return tProb

@cython.boundscheck(False)
@cython.wraparound(False)
def getInitProb(np.ndarray[np.float_t,ndim=1,negative_indices=False] I):
    cdef np.ndarray[np.float_t,ndim=1,negative_indices=False] ip = np.zeros((3),dtype=np.float)
    cdef float sum_k = I[0]+I[1]+I[2]
    cdef int i
    cdef float epsilon = 0.001
    if sum_k != 0.0:
       for i in range(0,3):
           ip[i] = I[i]/sum_k
    else:
       ip[0], ip[1], ip[2] = 1.0-2.0*epsilon,epsilon,epsilon
    return ip

@cython.boundscheck(False)
@cython.wraparound(False)
def new_eProb(np.ndarray[np.float_t,ndim=2,negative_indices=False] emm):
    cdef np.ndarray[np.float_t,ndim=2,negative_indices=False] eProb_bw = np.zeros((4,4),dtype=np.float)
    cdef np.ndarray[np.float_t,ndim=1,negative_indices=False] sum_k = np.zeros((4),dtype=np.float)
    cdef int i = 0
    cdef int j = 0
    for i in range(0,4):
        for j in range(0,4):
            sum_k[i] += emm[i,j]
    for i in range(0,4):
        if sum_k[i] != 0.0:
           for j in range(0,4):
               eProb_bw[i,j] = emm[i,j]/sum_k[i]
        else:
           if i == 0:#A
              eProb_bw[i,0],eProb_bw[i,1],eProb_bw[i,2],eProb_bw[i,3] = 0.5,0.1667,0.1667,0.5-2.0*(0.1667)
           if i == 1:#C
              eProb_bw[i,0],eProb_bw[i,1],eProb_bw[i,2],eProb_bw[i,3] = 0.1667,0.5,0.1667,0.5-2.0*(0.1667)
           if i == 2:#T
              eProb_bw[i,0],eProb_bw[i,1],eProb_bw[i,2],eProb_bw[i,3] = 0.1667,0.1667,0.5,0.5-2.0*(0.1667)
           if i == 3:#G
              eProb_bw[i,0],eProb_bw[i,1],eProb_bw[i,2],eProb_bw[i,3] = 0.1667,0.1667,0.5-2.0*(0.1667),0.5
    return eProb_bw
