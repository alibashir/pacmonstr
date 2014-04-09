#!/usr/bin/env python
import os
import sys
#sys.path.insert(0, '/projects/HuPac/repartition_MAR29/PacmonSTR/pHmm')
#sys.path.insert(0, '/projects/HuPac/repartition_MAR29/PacmonSTR/psCount')
#sys.path.insert(0, './pHmm')
#sys.path.insert(0, './psCount')
#sys.path.insert(0, './pacModel')
#sys.path.insert(0, './pacIO')
import numpy as np
from pacIO.ReadMatcherIO import parseRm5
from pacModel.AlignmentHit import AlignmentHit
from pacIO.binFileIO import parsebinf
from pacModel.binLine import binLine
from pacModel.Sequence import revcomp as rc
#from ReadMatcherIO import parseRm5
#from AlignmentHit import AlignmentHit
#from binFileIO import parsebinf
#from binLine import binLine
#from Sequence import revcomp as rc
import time
import multiprocessing as mp
import dpFuncs_sw2 as dpF
from psCount.mleHmm import countsOnList as cL
from psCount.mleHmm import getGlobalProb as load
from pHmm.pHmmFwdLog2_cdfLB import pairHmmRun 


def dictforbindata(binData,binfile):
    """
    The objective of this function is to form a dictionary from a file with query id and corresponding TR event information.

    Parameters
    ----------
    binData = dictionary
              {<query Ids>: <[list of TR events that this query Id spans]>}

    binfile = string
              Name of the input file with TR events for each query per line

    Returns
    -------
    Void

    """

    for line in parsebinf(binfile):
        binData.setdefault(line.query_id,[]).append(line)
    print "Length of binData: ", len(binData)

def dictforqueries(queryAligns,m5file):
    """
    The objective of this function is to form a dictionary from a file with query id and corresponding reference alignment information 

    Parameters
    ----------
    queryAligns = dictionary
              {<query Ids>: <reference alignment information>}

    m5file = string
              Name of the input .m5 file with reference alignment information 

    Returns
    -------
    Void

    """

    for hit in parseRm5(m5file):
        queryAligns[hit.query_id] = hit

def tb_prefix(query,prefix,P,P_p,iPre):
    """
    The objective of this function is to traceback the Prefix score matrix and return alignement sequences for prefix in reference and query
 
    Paramters
    ---------
    query = array of int, shape(n_bases_q,)
            The query sequence information

    prefix = array of int, shape(n_bases_p,)
             The prefix sequence information for the reference

    P = array of int, shape(n_bases_q,n_bases_p)
        The prefix score matrix

    P-p = array of int, shape(n_bases_q,n_bases_p)
          The prefix pointer matrix

    iPre = int
           row index where traceback jumps from TR repeat matrix R into prefix score matrix P

    Returns
    -------
    q_pre_Aligned = list
                    A list ['prefix sequence in ref','prefix sequence in query'] 

    """
    pre = False
    j = np.size(prefix)
    i = iPre
    qSeq = map(chr,query)
    preSeq = map(chr,prefix)
    qList = []
    preList= []
    jLessThanZero = False
    while(j > 0):
      if j > 0: 
         if (P[i][j] > 0):
            if (P_p[i][j] == 0):
               qList.append(qSeq[i-1])
               preList.append(preSeq[j-1])
               i += -1
               j += -1
            if (P_p[i][j] == 1):
               qList.append(qSeq[i-1])
               preList.append("-")
               i += -1
            if (P_p[i][j] == 2):
               qList.append("-")
               preList.append(preSeq[j-1])
               j += -1
         else:
            jLessThanZero = True
            break
      else:
         jLessThanZero = True
         break
    if jLessThanZero == False:
       preList.reverse()
       qList.reverse()
       q_pre_Aligned = []
       q_pre_Aligned.append(''.join(preList))
       q_pre_Aligned.append(''.join(qList))
    else:
       q_pre_Aligned = ['','']
    return q_pre_Aligned

def tb_suffix(query,suffix,S,S_p):
    """
    The objective of this function is to traceback the Suffix score matrix.
 
    Paramters
    ---------
    query = array of int, shape(n_bases_q,)
            The query sequence information

    suffix = array of int, shape(n_bases_s,)
             The suffix sequence information for the reference

    S = array of int, shape(n_bases_q,n_bases_s)
        The suffix score matrix

    S_p = array of int, shape(n_bases_q,n_bases_s)
          The suffix pointer matrix

    Returns
    -------
    index_i_tR = int
                 The index in query where suffix traceback enters the repeat matrix R

    bestscoresuffix = int
                      The maximum smith-waterman score

    q_suf_Aligned = list
                    A list ['suffix sequence in ref','suffix sequence in query'] 

    """

    temp = 0
    i = 0
    j = np.size(suffix)
    nq = np.size(query)
    qSeq = map(chr,query)
    sufSeq = map(chr,suffix)
    qList = []
    sufList= []
    for k in range(nq+1): # Finds the maximal value for the suffix to begin tb
        if temp < S[k][j]:
           temp = S[k][j]
           i = k
    bestscoresuffix = temp
    pre = False
    index_i_tR = 0
    jLessThanZero = False
    while (not(pre)):
      if j > 0:
         if (S[i][j] > 0):
            if (S_p[i][j] == 1):
               qList.append(qSeq[i-1])
               sufList.append("-")
               i += -1
            elif (S_p[i][j] == 2):
               qList.append("-")
               sufList.append(sufSeq[j-1])
               j += -1
            elif (S_p[i][j] == 0):
               qList.append(qSeq[i-1])
               sufList.append(sufSeq[j-1])
               i += -1
               j += -1
            elif (S_p[i][j] == 3):
               qList.append(qSeq[i-1])
               sufList.append(sufSeq[j-1])
               j += -1
               i += -1
               pre = True
            index_i_tR = i
         else :
            index_i_tR = 0
            jLessThanZero = True
            break
      else:
         index_i_tR = 0
         jLessThanZero = True
         break
    if jLessThanZero == False:
       sufList.reverse()
       qList.reverse()
       q_suf_Aligned = []
       q_suf_Aligned.append(''.join(sufList))
       q_suf_Aligned.append(''.join(qList))
       #print "q_suf_Aligned:\n",q_suf_Aligned[0],"\n",q_suf_Aligned[1] 
    else:
       q_suf_Aligned = ['','']
    return index_i_tR, bestscoresuffix,q_suf_Aligned

def finalTrace(qArray,tR,T,T_p,tR_rowMax,index_i_tR):
    """
    The objective of this function is to traceback the repeat score matrix and identify the TR sequence in query
 
    Paramters
    ---------
    qArray = array of int, shape(n_bases_q,)
            The query sequence information

    tR = string 
         The TR element  sequence information for the reference

    T = array of int, shape(n_bases_q + 1,n_bases_tR + 1)
        The repeat score matrix

    T_p = array of int, shape(n_bases_q + 1,n_bases_tR + 1)
          The repeat pointer matrix

    tR_rowMax = numpy array of int, shape(n_bases_q + 1, 2)
                stores the information on the maximum smith-waterman score for each row

    index_i_tR = int
                 The index from where the traceback starts in the repeat score matrix, T

    Returns
    -------
    number_tR = float
                Estimated number of TR elements in the query based on dynamic programming output

    repeatAligned = list
                    A list ['tandem repeat sequence in ref','tandem repeat sequence in query'] 

    scoretRentry = int
                   Score when entered into the repeat matrix, R

    scoretRexit = int
                  Score when exited from the repeat matrix into prefix score matrix

    naiveTR = float
              naive estimation of number of repeat elements based on estimated prefix and suffix boundaries

    """

    i = index_i_tR #This is from where the traceback starts
    #print "i index in TR: ", i
    if (i > 0):
      j = tR_rowMax[i][1]
      scoretRentry = T[i][j]
      t = len(tR)
      pre = False
      trScoreNeg = False
      query = map(chr,qArray)
      qList = []
      tList = []
      while (not (pre)):
        if j > 0:
          if (T[i][j] > 0):
            if (T_p[i][j] == 0):
               qList.append(query[i-1])
               tList.append("-")
               i += - 1
            elif (T_p[i][j] == 1):
               qList.append("-")
               tList.append(tR[j-1])
               j += - 1
            elif (T_p[i][j] == 2):
               qList.append(query[i-1])
               tList.append(tR[j-1])
               i += - 1
               j += - 1
            elif (T_p[i][j] == 3):
               qList.append(query[i-1])
               tList.append(tR[j-1])
               for k in range(0,j-1): #open up gaps on query
                   qList.append("-")
                   tList.append(tR[k])
               i += - 1
               j = t
            elif (T_p[i][j] == 4):
               qList.append(query[i-1])
               tList.append(tR[j-1])
               i += - 1
               pre = True
               trSinQ = i
               scoretRexit = T[i+1][j]
          else:
              trScoreNeg = True
              print "Tandem repeat score <= 0"
              break
        else:
            trScoreNeg = True
            break
      #print "trScoreNeg: ",trScoreNeg
      if trScoreNeg == False: #Repeat matrix score was never less than zero, implies all is well
         #print "tandem repeat score > 0"
         tList.reverse()
         qList.reverse()
         repeatAligned = []
         repeatAligned.append(''.join(tList))
         repeatAligned.append(''.join(qList))

         tempTList = []
         for i in range(0, len(tList)):
            if tList[i] == 'A' or tList[i] == 'C' or tList[i] == 'G' or tList[i] == 'T':
               tempTList.append(tList[i])

         tRLength = len(tempTList)
         str_tR = "".join(tR)
         str_tempTList = "".join(tempTList)
         if (t < tRLength):
              if t >= 1:
                 number_tR = (tRLength)/float(t) #This is the estimated number of tanderm repeat elements from the dynamic programming method
                 naiveTR = index_i_tR - trSinQ #naive number of tandem repeat elements based on boundary estimation
              else:
                 number_tR = -30000.0
                 naiveTR = -30000.0
                 repeatAligned = ['','']
                 scoretRentry = 0
                 scoretRexit = 0
         else:
           number_tR = -20000.0
           naiveTR = -20000.0
           repeatAligned = ['','']
           scoretRentry = 0
           scoretRexit = 0
      else:
          number_tR = -40000.0
          naiveTR = -40000.0
          repeatAligned = ['','']
          scoretRentry = 0
          scoretRexit = 0
    else:
      number_tR = -10000.0
      naiveTR = -10000.0
      repeatAligned = ['','']
      scoretRentry = 0
      scoretRexit = 0
    return number_tR, repeatAligned, scoretRentry, scoretRexit, naiveTR
 
def alignRegions(query,pre_suf_tR):
    """
    The objective of this function is to identify the TR boundaries within the reads (Section 2.1)

    Parameters
    ---------- 
    query = string 
           The read sequence

    pre_suf_tR = List
                 A list of strings: [<prefix_seq>, <suf_seq>, <tr_seq>]

    Returns
    -------  
    #return number_tR, repeatAligned, Prefix_BestScore, Suffix_BestScore, tR_maxScore, tR_exitScore, alignSufxQ, alignPrfxQ, naiveTR
    number_tR = float
                number of tandem repeats in query as estimated by the 3-stage dynamic programming algorithm

    repeatAligned = list
                    A list of strings: [<TR align sequence in ref>, <TR align sequence in query>]

    Prefix_BestScore = int
                       Smith-waterman score of prefix score matrix, P (Prefix)

    Suffix_BestScore = int
                       Smith-waterman score in suffix score matrix, S (Prefix + TR + Suffix)

    tR_maxScore = int
                  Smith-waterman score in TR matrix, R (Prefix + TR)

    tR_exitScore = int
                   Smith-waterman score in TR matrix where alignment jumps to prefix matrix, P

    alignSufxQ = list
                 A list of Aligned suffix sequence in read ['suffix align sequence in reference','suffix align sequence in query']

    alignPrfxQ = list
                 A list of aligned prefix sequence in read ['prefix align sequence in reference','prefix align sequence in query']

    naiveTR = float
              number of tandem repeats in query as estimated by naive estimation of prefix and suffix boundaries in read
 
    Notes
    -----
    This function implements section 2.1 from the paper. For speed ups the main dynamic programming parts have been implemented in cython.
    """
    temp_trfSeq = []
    prefix = []
    prefix = list(pre_suf_tR[0])
    len_prefix = len(prefix)
    pArray = np.arange(len_prefix,dtype=np.int_)
    tempP = map(ord,prefix)
    tempP2 = np.asarray(tempP)
    pArray = tempP2 #A numpy array of the PREFIX
#    print "pArray: ", pArray

    suffix = []
    suffix = list(pre_suf_tR[1])
    len_suffix = len(suffix)
    sArray = np.arange(len_suffix,dtype=np.int_)
    tempS = map(ord,suffix)
    tempS2 = np.asarray(tempS)
    sArray = tempS2 #A numpy array of the SUFFIX
#    print "sArray: ", sArray

    tR = []
    tR = list(pre_suf_tR[2])
    len_tR = len(tR)
    tRArray = np.arange(len_tR,dtype=np.int_)
    temptR = map(ord,tR)
    temptR2 = np.asarray(temptR)
    tRArray = temptR2 #A numpy array of the TR element
#    print "tRArray: ", tRArray

    #initialization for prefix
    preScore = np.zeros((len(query)+1,len_prefix+1),dtype=np.int_)
    # Added in initialization for prefix
    preMaxScore = np.zeros((len(query)+1,len_prefix+1),dtype=np.int_)
    preMaxScore.fill(-1)
    for j in xrange (1, len_prefix+1):
        preScore[0][j] = -5*j

    Prefix_BestScore = 0
    Prefix_BestScoreIndex = 0
    Prefix_BestScore, Prefix_BestScoreIndex = dpF.sw_prefix(query,pArray,preScore,preMaxScore) #Calling sw_prefix function (implemented in Cython)

    #initialization for TR
    tRScore = np.zeros((len(query)+1,len_tR+1),dtype=np.int_)
    for j in xrange(1,len_tR+1):
        tRScore[0][j] = -10000000

    tRMaxScore = np.zeros((len(query)+1,len_tR+1),dtype=np.int_)-1
    tR_rowMax = np.zeros((len(query)+1,2),dtype=np.int_)
    dpF.sw_tR_simple(query,tRArray,tRScore,tRMaxScore,len_prefix,preScore,tR_rowMax) #Calling sw_tR_simple function (implemented in Cython)

    #initialization for suffix
    index_i_tR = 0
    Suffix_BestScore = 0
    sufScore = np.zeros((len(query)+1,len_suffix+1),dtype=np.int_)
    for j in xrange(1,len_suffix+1):
        sufScore[0][j] = -10000000

    sufMaxScore = np.zeros((len(query)+1,len_suffix+1),dtype=np.int_)-1
    dpF.sw_suffix(query,sArray,sufScore,sufMaxScore,tR_rowMax) #Calling sw_suffix function (implemented in Cython)
    index_i_tR,Suffix_BestScore,alignSufxQ = tb_suffix(query,sArray,sufScore,sufMaxScore)

    number_tR = 0.0
    naiveTR = 0.0
    tR_maxScore = 0
    repeatAligned = []
    number_tR, repeatAligned, tR_maxScore, tR_exitScore, naiveTR = finalTrace(query,tR,tRScore,tRMaxScore,tR_rowMax,index_i_tR) #Traceback from the suffix matrix, S and into TR repeat matrix, R 

    alignPrfxQ = tb_prefix(query,pArray,preScore,preMaxScore,Prefix_BestScoreIndex) #Here goes the prefix matrix traceback. here query is in int, but tR is in string representation
    print "exit finalTrace..."
    return number_tR, repeatAligned, Prefix_BestScore, Suffix_BestScore, tR_maxScore, tR_exitScore, alignSufxQ, alignPrfxQ, naiveTR

#####

def calculateRegions_tR(flag,target, rev_target, trf_start, trf_end, target_start, target_end,refTrf_seq, flankLen):
    """
    The objective of this function is to extract out the prefix, suffix & TR element sequences from the reference alignment and TR event information passed as input.

    Parameters
    ----------
    flag = bool type (True or False)
           Indicates positive or negative strand of the reference where the query aligned

    target = string
             Reference sequence information passed from an alignment object  
    
    rev_target = string
                 Reverse complement of the reference sequence information in the form of a string

    trf_start = int
                TR event start position or coordinates on the reference
   
    trf_end = int
              TR event end position or coordinates on the reference

    target_start = int
                   Position where query starts alignment with the reference sequence

    target_end = int
                 Position where query ends alignment with the reference sequence
 
    refTrf_seq = string
                 The TR element in the reference

    flankLen = int
               The length of flank region around the TR event

    Returns
    -------
    listPreSuf_tR = List
                    A list of three elements [prefix flank region; suffix flank region; reference TR element sequence]
    
    """
    listPreSuf_tR = []
    tempTarget = []
    prefix_tR = ""
    suffix_tR = ""

    lengthOfPreX = trf_start - target_start
    lengthOfSufX = target_end - target_start
    startSufX = trf_end - target_start
    len_target = len(target)
    len_revtarget = len(rev_target)

    if (len_target < lengthOfSufX):
       listPreSuf_tR.append('E')
    else:
       if flag == True:
          tempRevTarget = rev_target
          rev_prefix_tR = []
          for j in xrange(0,lengthOfPreX):
              rev_prefix_tR.append(tempRevTarget[j])
          suffix_tR = rc(''.join(rev_prefix_tR))
          rev_suffix_tR = []
          for k in xrange(startSufX, lengthOfSufX):
              rev_suffix_tR.append(tempRevTarget[k])
          prefix_tR = rc(''.join(rev_suffix_tR))

          refTrf_seq = rc(refTrf_seq)
       else:
          tempTarget = target
          prefix_tR = []
          for j in xrange(0,lengthOfPreX):
              prefix_tR.append(tempTarget[j])

          suffix_tR = []
          for k in xrange(startSufX, lengthOfSufX):
              suffix_tR.append(tempTarget[k])
          prefix_tR = ''.join(prefix_tR)
          suffix_tR = ''.join(suffix_tR)
       listPreSuf_tR.append(prefix_tR[-flankLen:])
       listPreSuf_tR.append(suffix_tR[0:flankLen])
       listPreSuf_tR.append(refTrf_seq)
    return listPreSuf_tR

def accAlign(sufAlign):
    sufMatch = 0
    sufL = len(sufAlign[0])
    for i in range(0,sufL):
        if sufAlign[0][i] == sufAlign[1][i]:
           sufMatch += 1
    sufAcc = sufMatch/float(sufL)
    return sufAcc

def getAccFromAlign(listOfAligned,repeatAligned):
    preAlign = listOfAligned[1]
    preL = len(preAlign[0])
    preMatch = 0
    sufAcc = accAlign(listOfAligned[0])
    preAcc = accAlign(listOfAligned[1])
    trAcc = accAlign(repeatAligned)
    return preAcc, sufAcc, trAcc

def align_worker(ListOfQs,m5Dict,binData,fn):
    """
    This function runs the entire pipeline for pacmonSTR.

    Parameters
    ----------
    This function takes in 4 inputs: 
    ListOfQs = List
               List of query Ids in the form of strings, [qId_1, qId_2,...,qId_k]
  
    m5Dict = Dictionary
             A dictionary with alignment information for each query_id (m5Dict = {<qId_i:Alignment_Object_i>})

    binData = Dictionary
              A dictionary with TR event information for each query_id (binData = {<qId_i:[TR_Event_Object_1, TR_Event_Object_1, ... , TR_Event_Object_i]})
              A query can span multiple TR events, and hence for every query we store the list of TR Event objects
    fn = string
         Name of the outfile where results would be written (fn = a string)

    Returns
    -------
    Void
    Writes the results of the pipeline to a file with the name 'fn'

    Notes
    -----
    The pipeline goes through the following key stages:
    A. Given an alignment between a query and the reference, it extracts the prefix and suffix strings from the reference.
    B. Identification of TR boundaries within reads Qi (Section 2.1).It forms the prefix (P), suffix (S) and repeat matrices (R) and goes through the 3-stage dynamic programming step to identify the tandem repeat sequence, q from Qi.
    C. Probabilistic TR resolution through pairHMM (Section 2.2 & 2.3). For each identified tandem repeat sequence, q,  within the read Qi, it estimates the expected TR multiplicity and the probable SV regions in q as a list of tuples, or, regionHigh .
    D. TR allele calling (Section 2.4): Given the list of estimated TRs for all the reads spanning a particular TR event, it uses GMM/AIC based model selection to call a homozygous or heterozygous event and runs POA to generates a consensus sequence for each cluster.
    """

    start1 = time.time()
    FinalList = [] #A list for storing output for each TR event spanned by the query, qId_i
    for query in ListOfQs: #Iterate through the list of queries
        hit = m5Dict[query]
        if hit.query_id in binData: #Check to see if that query has a corresponding TR event by checking to see if it is a key in the binData dictionary
           tempLine = binLine() #Instantiate a object of type binLine() (Basically this has the TR event information)
           loq = len(hit.QuerySeq)
           qSeq = ''.join(hit.QuerySeq)
           qArray = np.arange(loq,dtype=np.int_)
           tempQ = map(ord,hit.QuerySeq)
           tempQ2 = np.asarray(tempQ)
           qArray = tempQ2 #Initialized a numpy array of corresponding integers to the bases
           pre_suf_tR = [] #List which will store sequence information for the Prefix (P), Suffix (S) and TR repeat element
           flag = False # assume its in positive strand
           if hit.target_strand == 0:
              flag = True #1 means its the negative strand
           numInList = len(binData[hit.query_id])
           for k in range(0,numInList): #go through all the TR events for that query
               tempLine = binData[hit.query_id][k]
               if tempLine.prefixFlank >= 100 and tempLine.suffixFlank >= 100: #check to see if the flanks around the TR event satisfies a certain criterion (of 100bp here)
                   pre_suf_tR = calculateRegions_tR(flag,hit.TargetSeq,hit.RevTargetSeq,tempLine.refTrf_st,tempLine.refTrf_en,hit.target_start,hit.target_end,tempLine.refTrf_seq,flankLen) #Gets all the required sequences for prefix, suffix and TR element
                   if (len(pre_suf_tR)) == 1: #If this is size one then something is not right with this alignment. Report and exit
                      tempStr = "query has weird m5 alignment (check m5 file alignment)\n"
                      outList = [hit.query_id, tempStr]
                   else:
                      number_tR, repeatAligned, PreS, SufS, tR_S, tR_Ex, qSufAligned, qPreAligned, naiveTR = alignRegions(qArray, pre_suf_tR) #Runs the 3-stages dynamic programming step to identify the tandem repeat sequence, q
                      print "Done alignRegions"
                      if repeatAligned[1] != '': #Check if TR in q exists 
                         tList = repeatAligned[1] #This is the alignment string of the tandem repeat sequence in query
                         tempTList = []
                         for i in range(0, len(tList)): #here all the gaps or '-' are removed 
                             if tList[i] == 'A' or tList[i] == 'C' or tList[i] == 'G' or tList[i] == 'T':
                                tempTList.append(tList[i])
                         trInQseq = ''.join(tempTList) #This is the filtered string or the tandem repeat sequence in query, q 
                         preSeqInQ = ''.join([x if x != '-' else '' for x in qPreAligned[1]]) #all the '-' are removed for the identified prefix sequence in query
                         sufSeqInQ = ''.join([x if x != '-' else '' for x in qSufAligned[1]]) #all the '-' are removed for the identified suffix sequence in query
                         if len(preSeqInQ) >= 100 and len(sufSeqInQ) >= 100: #flanking prefix and suffix length are more than 100. Perform both Local and Global pHmm calculations
                            preFlank = preSeqInQ[-95:]
                            sufFlank = sufSeqInQ[0:95]
                            listOfAligned = [qSufAligned, qPreAligned]
                            preAcc, sufAcc, trAcc = getAccFromAlign(listOfAligned,repeatAligned) #get the alignment accuracy scores for the prefix and suffix regions. Needs to fulfill a particular cut-off for alignment, else that flank is not trustworthy
                            print "preAcc, sufAcc, trAcc ", preAcc, " ", sufAcc, " ", trAcc 
                            if preAcc >= 0.75 and sufAcc >= 0.75: # accuracy above 75% so that the flanks around the tR event are reliable 
                               tProb, eProb_M, eProb_X = cL(listOfAligned) #Local Parameters for pHmm (lambda_local parameter set). These are calculated for the local sequences. 
                               tProb1, eProb_M1, eProb_X1 = load(hit.target_id) #Global Parameters for pHmm (lambda_global parameter set). These are pre-computed and loaded here
                               print "GLocal pHmm start..."
                               numTRL, llrL, ProbMaxL, SamplesL, lowerBoundL, listOfDevL = pairHmmRun(trInQseq,pre_suf_tR[2],tProb,eProb_M,eProb_X) #This is the Local pairHMM run using the lambda_local parameter set
                               numTRG, llrG, ProbMaxG, SamplesG, lowerBoundG, listOfDevG = pairHmmRun(trInQseq,pre_suf_tR[2],tProb1,eProb_M1,eProb_X1) #This is the global pairHMM run using the lambda_global parameter set
                               print "Global psCount pHmm End"
                               print "number_tR: ",number_tR," numTRL: ",numTRL," numTRG: ",numTRG
                               if ProbMaxL > ProbMaxG: #Check which run of pHmm gives higher probability of P(O|lambda). Select those parameters respectively 
                                  numTR,llr,ProbMax,Samples,lowerBound,listOfDev = numTRL, llrL, ProbMaxL, SamplesL, lowerBoundL, listOfDevL
                               else:
                                  numTR,llr,ProbMax,Samples,lowerBound,listOfDev = numTRG, llrG, ProbMaxG, SamplesG, lowerBoundG, listOfDevG
                            else:
                               numTR, llr, ProbMax, Samples, lowerBound, listOfDev = -1,0,-1.0*float('inf'),[],0,[] 
                         else: #Don't want q which don't span the TR event greater than 100 base pairs
                            preFlank, sufFlank = '', ''
                            numTR, llr, ProbMax, Samples, lowerBound, listOfDev = -1,0,-1.0*float('inf'),[],0,[]
                            preAcc,sufAcc,trAcc = 0.0,0.0,0.0
                      else: #TR is q doesn't exists...pass default values
                         numTR, llr, ProbMax, Samples, lowerBound, listOfDev = -1,0,-1.0*float('inf'),[],0,[]
                         repeatAligned = ["Error","Error"]
                         preFlank, sufFlank = '', '' #We don't want pre and suf in q sequences if there is no TR in q detected
                         preAcc,sufAcc,trAcc = 0.0,0.0,0.0
                      outList = [hit.query_id,hit.target_id,hit.target_start,hit.target_end,pre_suf_tR[2],number_tR,tempLine.refTrf_copyNum,preAcc,sufAcc,trAcc,tR_Ex-tR_S,repeatAligned[0],repeatAligned[1],tempLine.refTrf_st,tempLine.refTrf_en,tempLine.refTrf_repeat,naiveTR/float(tempLine.refTrf_repeat),numTR,llr,ProbMax,Samples,lowerBound,listOfDev,preFlank,sufFlank]
                   FinalList.append(outList)
               else:
                   tempStr = "query flanks around TR event less than 100 bp\n"
                   outList = [hit.query_id, tempStr]
                   FinalList.append(outList)
        else:
           tempStr = "query does not exist in binned file\n"
           outList = [hit.query_id, tempStr]
           FinalList.append(outList)
    end1 = time.time()
    print "Duration within worker: ",end1-start1
    writeFile(FinalList,fn) #Write the final result to the file, 'fn' in the input run directory
    print "worker ends for query" 

def writeFile(x,fn):
    """
    The objective of this function is to write the results to a file name: 'fn'

    Parameters
    ----------
    This function takes in two inputs:
    x = A list of list with all the output results
    fn = name of the file to write the output result
   
    Returns
    -------
    Void

    Notes: Iterates through the list x and converts each list into a string and writes to a file
    """
    print "printing fine..."
    outfile = open(fn,'w') 
    for result in x:
        strResult = ' '.join(map(str,result))
        outfile.write("%s\n"%(strResult))
        outfile.flush()
    outfile.close()


if __name__ == '__main__': 
    #'''USAGE: python pacMonStr_V0.py <>.m5 <>.binned 1000 <numCores> <outputDirectory>'''
    start = time.time()
    inm5fn = sys.argv[1] #input 1: Alignment file in m5 format from blasr
    tempstr = inm5fn.rstrip("\n").split("/")
    ver = '.PacV1'

    outdir = sys.argv[5] #input 5: Output directory where the results would be written

    outfn = tempstr[-1].split(".")[0] + ver
    outStr = outdir + "/" + outfn + ".out"

    inbinfn = sys.argv[2] #input 2: .binned_anchors file which has reads filtered based on TR events with all information on the TR event

    flankLen = int(sys.argv[3]) #input 3: length of region flanking TR (epsilon in the paper)

    ncores = int(sys.argv[4]) #input 4:  number of cpus to use

    binData = {} #binData is a global dictionary where each query_id is the key and the related TR event it spans is the value
    dictforbindata(binData,inbinfn) #<query_id:filtered queries for a TR event with its relevant information> entries
    print "No of entries in binData dict: ",len(binData)

    ListOfQueryAligns = []
    queryAligns = {} #queryAligns is a dictionary where each query_id is the key and the related alignment with a reference (from .m5 file) is the value
    dictforqueries(queryAligns,inm5fn)
    NumOfQ = len(queryAligns)
    print "No of entries in queryAligns dict: ", NumOfQ

    
    num = int(NumOfQ/(NumOfQ*0.05)) # 5% of number of queries, so around 100
    if NumOfQ > num:
       chunkDict = int(NumOfQ/num) #a parameter which is used to partition the dictionary, queryAligns, and passes that partition to a processor via pool
    else:
       chunkDict = 1
    ListOfQueries = queryAligns.keys()
    chunkedListQueries = [ListOfQueries[i::chunkDict] for i in range(chunkDict)] #make the partitions to be sent to each pool 
    pool = mp.Pool(ncores)
    for i in range(0,len(chunkedListQueries)): #send each partition via pool for processing
         fn = outStr +"_" +str(i+1)
         print "fn: ",fn
         pool.apply_async(align_worker,(chunkedListQueries[i],queryAligns,binData,fn)) #align_worker takes each chunk and processes it before exiting. Eveything is done in this function.
    pool.close() #Close all the pool objects and join them. Waits for all the processes to be finished before exiting
    pool.join()
    print "Done...exiting PacmonSTR"
