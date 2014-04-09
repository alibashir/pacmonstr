from model.Sequence import revcomp as rc
import sys

class AlignmentHit:
    def __len__( self ):
        return self.query_end - self.query_start
    
    def __init__( self ):
        self.query_length, self.query_id, self.query_start, \
            self.query_end, self.query_strand = 0, None, 0, 0, None
        self.target_length, self.target_id, self.target_start, \
            self.target_end,self.target_strand = 0, None, 0, 0, None
        self.score = 0
        self.aligned = ""
        self.alignedQuery = ""
        self.alignedTarget = ""
        self.QuerySeq = []
        self.QueryStr = ""
        self.QueryStrOrg = ""
        self.TargetSeq = []
        self.TargetStr = ""   
        self.TargetStrOrg = ""   
        self.RevQuerySeq = []
        self.RevQueryStr = ""
        self.RevTargetSeq = []
        self.RevTargetStr = ""  
 
    def __str__( self ):
        return 'agar: %d %d %s %d %d %s %s %d %d %s %d' % \
            ( self.query_length, self.target_length, \
              self.query_id, self.query_start, self.query_end, \
              self.query_strand, self.target_id, self.target_start, \
              self.target_end, self.target_strand, self.score )

    def getAligned(self):
        return self.QueryStrOrg,self.TargetStrOrg

    def revcomp( self ):
        self.aligned = self.aligned[::-1]
        self.alignedQuery = rc(self.alignedQuery)
        self.alignedTarget = rc(self.alignedTarget)
        if self.target_strand == 1:
            self.target_strand = 0
        else:
            self.target_strand = 1
    
    def targetPosToIndex (self, pos):
        curr_pos = self.target_start
        target = self.alignedTarget
        if self.target_strand == 1:
            target = target[::-1]
        index = 0
        while curr_pos < pos:
            if target[index] != "-":
                curr_pos += 1
            index += 1
        return index


    def targetPosToIndices (self, spos, epos):
        curr_pos = self.target_start
        target = self.alignedTarget
        if self.target_strand == 1:
            target = target[::-1]
        index = 0
        sindex = 0
        eindex = 0
        while curr_pos < epos:
            if curr_pos == spos:
                sindex = index
            if target[index] != "-":
                curr_pos += 1
            index += 1
        return sindex, index


    def getUpstreamAnchor (self, x, anchor):
        index_e = x
        index_s = x-anchor
        errors = self.aligned[index_s:index_e].count("*")
        enderror = 0
        while index_s >= 0:
            if errors == 0:
                return index_s
            if self.aligned[index_e-1] == "*":
                errors += -1
            index_s += -1
            index_e += -1
            if self.aligned[index_s] == "*":
                errors += 1

        return -1


    def getDownstreamAnchor (self, x, anchor):
        index_s = x
        index_e = x+anchor
        errors = self.aligned[index_s:index_e].count("*")
        enderror = 0
        while index_e < len(self.aligned):
            if errors == 0:
                return index_e
            if self.aligned[index_s] == "*":
                errors += -1
            index_s += 1
            index_e += 1
            if self.aligned[index_e-1] == "*":
                errors += 1

        return -1
            
        

    def getBetweenHistDist (self, hit):
        qdist, tdist = -sys.maxint-1, -sys.maxint-1
        if hit.target_id != self.target_id:
            #return -sys.maxint-1, -sys.maxint-1, -sys.maxint-1, -sys.maxint-1
            return qdist, tdist
        #qs_min, qs_max, ts_min, ts_max = -1, -1, -1, -1
        if (hit.query_start < self.query_start and 
            hit.query_end < self.query_end):
            qdist = self.query_start-hit.query_end
            tdist = self.target_start-hit.target_end
        elif (self.query_start < hit.query_start and 
              self.query_end < hit.query_end):
            qdist = hit.query_start-self.query_end
            tdist = hit.target_start-self.target_end
            
        return qdist, tdist

        
        
        
        
