import numpy as np
cimport numpy as np

np.import_array()

def sw_prefix(np.ndarray[np.int_t,ndim=1] query, np.ndarray[np.int_t,ndim=1] prefix, np.ndarray[np.int_t,ndim=2] preScore, np.ndarray[np.int_t,ndim=2] preMaxScore):
    #declarations:
    cdef int nq, npre, pre, m

    nq = np.size(query)
    npre = np.size(prefix)
    pen = -5

    cdef np.ndarray[np.int_t, ndim=1] vals = np.zeros((3),dtype=np.int)

    cdef int bestscr  
    cdef int bestscrIndex  
    cdef int col

    #Py_ssize_t is a signed integer type
    cdef Py_ssize_t i, j, k
    cdef int temp, argtemp
    for i in range(1,nq+1):
        for j in range(1,npre+1):
            m = -6
            if query[i-1] == prefix[j-1]:
               m = 5
            vals[0] = preScore[i-1,j-1] + m
            vals[1] = preScore[i-1,j] + pen
            vals[2] = preScore[i,j-1] + pen
            #preScore[i,j] = np.max(vals)
            #preMaxScore[i,j] = np.argmax(vals)
            temp = vals[0]
            argtemp = 0
            for k in [1,2]:
                if vals[k] > temp:
                   temp = vals[k]
                   argtemp = k
            preScore[i,j] = temp
            preMaxScore[i,j] = argtemp
    bestscr = 0
    col = npre
    for i in range(nq+1):
        if bestscr < preScore[i,col]:
           bestscr = preScore[i,col]
           bestscrIndex = i
    return bestscr, bestscrIndex    


def sw_tR_simple(np.ndarray[np.int_t,ndim=1] query,np.ndarray[np.int_t,ndim=1] tR, np.ndarray[np.int_t,ndim=2] T, np.ndarray[np.int_t,ndim=2] T_p, int len_prefix, np.ndarray[np.int_t,ndim=2] P, np.ndarray[np.int_t,ndim=2] T_rowMax):
    cdef int p, t, nq, pen, m
    p = len_prefix
    t = np.size(tR)
    nq = np.size(query)
    pen = -5

    cdef np.ndarray[np.int_t, ndim=1] vals = np.zeros((5),dtype=np.int)

    #Py_ssize_t is a signed integer type
    cdef Py_ssize_t i, j, k
    cdef int temp, argtemp
    for i in range(1,nq+1):
        for j in range(1,t+1):
            m = -6
            if query[i-1] == tR[j-1]:
                m = 5
            vals[0] = T[i-1,j] + pen
            vals[1] = T[i,j-1] + pen
            vals[2] = T[i-1,j-1] + m
            vals[3] = T[i-1,t] + m + (j-1)*pen
            vals[4] = P[i-1,p] + m
            #T[i,j] = np.amax(vals)
            #T_p[i,j] = np.argmax(vals)
            temp = vals[0]
            argtemp = 0
            for k in [1,2,3,4]:
                if vals[k] > temp:
                   temp = vals[k]
                   argtemp = k
            T[i,j] = temp
            T_p[i,j] = argtemp

    cdef int temp0
    cdef int temp_tR
    #Py_ssize_t is a signed integer type
    cdef Py_ssize_t k0, l
    for k0 in range(nq+1):
        temp0 = T[k0,0]
        temp_tR = 0
        for l in range(1,t+1):
            if temp0 < T[k0,l]:
               temp0 = T[k0,l]
               temp_tR = l
        T_rowMax[k0,0] = temp0
        T_rowMax[k0,1] = temp_tR


def sw_suffix(np.ndarray[np.int_t,ndim=1] query,np.ndarray[np.int_t,ndim=1] suffix,np.ndarray[np.int_t,ndim=2] sufScore, np.ndarray[np.int_t,ndim=2] sufMaxScore, np.ndarray[np.int_t,ndim=2] tR_rowMax):
    cdef int nq, nsuf, pen, m
    nq = np.size(query)
    nsuf = np.size(suffix)
    pen = -5

    cdef np.ndarray[np.int_t, ndim=1] vals = np.zeros((3),dtype=np.int)

    #Py_ssize_t is a signed integer type
    cdef Py_ssize_t i, j, k
    cdef int temp, argtemp
    for i in range(1,nq+1):
        m = -6
        if query[i-1] == suffix[0]:
           m = 5
        sufScore[i-1,1] = tR_rowMax[i-1,0] + m
        sufMaxScore[i-1,1] = 3
        for j in range(2,nsuf+1):
            m = -6
            if query[i-1] == suffix[j-1]:
               m = 5
            vals[0] = sufScore[i-1,j-1] + m
            vals[1] = sufScore[i-1,j] + pen
            vals[2] = sufScore[i,j-1] + pen
            #sufScore[i,j] = np.max(vals)
            #sufMaxScore[i,j] = np.argmax(vals)
            temp = vals[0]
            argtemp = 0
            for k in [1,2]:
                if vals[k] > temp:
                   temp = vals[k]
                   argtemp = k
            sufScore[i,j] = temp
            sufMaxScore[i,j] = argtemp
