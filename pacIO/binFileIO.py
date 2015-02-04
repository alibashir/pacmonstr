from itertools import *
from pacModel.binLine import binLine

def parsebinf(file):
    for line in open(file):
        values = line.rstrip("\n").split(" ")
        hit = binLine()
        hit.query_id, hit.target_start, hit.target_end, hit.refTrf_st, hit.refTrf_en = values[0], int(values[1]), int(values[2]), int(values[3]),int(values[4])
        hit.refTrf_repeat, hit.refTrf_copyNum, hit.refTrf_seq = float(values[5]), float(values[6]), values[7]
        hit.prefixFlank, hit.suffixFlank = int(values[8]), int(values[9])
        yield hit
    

