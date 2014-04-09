from itertools import *
from model.AlignmentHit import AlignmentHit

def parseRm5( file ):
    """Parses readmatcher -printFormat 5 output into AlignmentHit objects"""

    for line in open(file):
        #print line
        values = line.rstrip("\n").split(" ")
        hit = AlignmentHit()
        hit.query_id, hit.query_length, hit.query_start, hit.query_end, hit.query_strand = values[0], int(values[1]), int(values[2]), int(values[3]), values[4]

        hit.target_id, hit.target_length, hit.target_start, hit.target_end, hit.target_strand = values[6], int(values[7]), int(values[8]), int(values[9]), values[10]
        hit.score = -1*int(values[11])
        #target_id = values[6] because there is white space in the m5 file before target_id
 
        hit.query_id = "/".join(hit.query_id.split("/")[0:3])
 
        hit.alignedQuery = values[17]
        hit.QueryStrOrg = hit.alignedQuery
        tempQList = []
        tempQList = list(hit.alignedQuery)
        for i in range(0, len(tempQList)):
            if tempQList[i] == 'A' or tempQList[i] == 'C' or tempQList[i] == 'G' or tempQList[i] == 'T':
               hit.QuerySeq.append(tempQList[i])
        hit.QueryStr = ''.join(hit.QuerySeq)       
        
        hit.alignedTarget = values[19]
        hit.TargetStrOrg = hit.alignedTarget
        tempTList = []  
        tempTList = list(hit.alignedTarget)
        for i in range(0, len(tempTList)):
            if tempTList[i] == 'A' or tempTList[i] == 'C' or tempTList[i] == 'G' or tempTList[i] == 'T':
               hit.TargetSeq.append(tempTList[i])
        hit.TargetStr = ''.join(hit.TargetSeq)       
        
        hit.aligned = values[18]      
        hit.line = line
        
        #print hit.target_strand
        if hit.target_strand == "+":
            hit.target_strand = 0
        else: 
            hit.target_strand = 1

        if hit.query_strand == "+":
            hit.query_strand = 0
        else: 
            hit.query_strand = 1
           
        hit.revcomp()

        tempRevQList = []
        tempRevQList = list(hit.alignedQuery)
        for i in range(0, len(tempRevQList)):
            if tempRevQList[i] == 'A' or tempRevQList[i] == 'C' or tempRevQList[i] == 'G' or tempRevQList[i] == 'T':
               hit.RevQuerySeq.append(tempRevQList[i])
        hit.RevQueryStr = ''.join(hit.RevQuerySeq)

        tempRevTList = []
        tempRevTList = list(hit.alignedTarget)
        for i in range(0, len(tempRevTList)):
            if tempRevTList[i] == 'A' or tempRevTList[i] == 'C' or tempRevTList[i] == 'G' or tempRevTList[i] == 'T':
               hit.RevTargetSeq.append(tempRevTList[i])
        hit.RevTargetStr = ''.join(hit.RevTargetSeq)

        yield hit

def parseRm4( file ):
    """Parses readmatcher -printFormat 4 output into AlignmentHit objects"""

    for line in open(file):
        fields = line.rstrip("\n").split(" ")
        hit = AlignmentHit()
        hit.query_id = fields[0]
        hit.target_id = fields[1]
        hit.score = - int(fields[2])
        hit.pctidentity = float(fields[3])
        hit.query_strand = "+" if fields[4] == "0" else "-"

        hit.query_start =  int(fields[5])
        hit.query_end  =   int(fields[6])
        hit.query_length = int(fields[7])

        # for negative strand readMatcher query coords are reported on reverse 
        # complement of the sequence. For the alignmenthit they need to be 
        # reported on the forward strand
        if hit.query_strand == "-":
            hit.query_end, hit.query_start = hit.query_length - hit.query_start, hit.query_length - hit.query_end 

        hit.target_strand = "+" if fields[8] == "0" else "-"
        hit.target_start  = int(fields[9])
        hit.target_end    = int(fields[10])
        hit.target_length = int(fields[11])
        
        if hit.target_strand == "-":
            hit.target_end, hit.target_start = hit.target_length - hit.target_start, hit.target_length - hit.target_end 
        
        #if hit.target_strand == hit.query_strand: # always report strand for query
        #    hit.target_strand = hit.query_strand = "+"
        #else:
        #    hit.query_strand = "-"
        #    hit.target_strand = "+"

        hit.alignedLength = abs(hit.target_end-hit.target_start)
        yield hit

