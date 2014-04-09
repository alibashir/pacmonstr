rcDict = dict(zip('ACGTWSMKRYBDHVNacgtwsmkrybdhvn-','TGCAWSKMYRVHDBNtgcawskmyrvhdbn-'))

def revcomp(seq):
   return ''.join(map(lambda x: rcDict[x],seq[::-1]))

