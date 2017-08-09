import os
import distance

from cutadapt.align import Aligner
from cutadapt.adapters import ANYWHERE


## extracts the barcodes and the UMI ##
## 2 step cutadapt to make sure that we find the linker positions - could be improved 
## could not get BWApy to work a API level wrapper around an aligner returning a CIGAR string would be much better
## 

PS = 5
PE = 48
L1A = "TAGCCATCGCATTGC"
L2A = "TACCTCTGAGCTGAA"
L2B = "TACCTCTTATCTCTT"
PRIME5 = "ACGGAC"

L1s = { "L1A" : L1A }
L2s = { "L2A" : L2A, "L2B" : L2B }


def getLinker1(name):
    return L1s[name]
    
def getLinker2(name):
    return L2s[name]

def getLLinker(name1, name2):
    return getLinker1(name1)+"NNNNNN"+getLinker2(name2)


def loadBarcodes():
    barcodes = []
    with open(os.path.dirname(os.path.abspath(__file__))+"/barcodes.txt") as b:
        for l in b:
            barcodes.append(l.strip())
    return barcodes


BARCODES = loadBarcodes()
BARCODESMAP = dict(zip(BARCODES, BARCODES))
EXTRACTCOUNTER = {"noLL": 0, "noL1orL2": 0,
                "L1not6fromL2": 0, "noACGorGAC": 0, "extracted": 0}
BCCOUNTER = {"exact": 0, "mm1": 0, "none": 0}
PRIME5COUNTER = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}

# mutation rate full linker: length is 6 + 15+6+15 + 6 = 48: 0.05
# mutation rate linker part: 0.1


def findLinker(readPart, linker, mut):
    aligner = Aligner(readPart, mut, ANYWHERE,
                      wildcard_ref=True, wildcard_query=True)
    l = aligner.locate(linker)
    if l is None:
        return None
    rStart = l[0]
    rEnd = l[1]
    matches = l[4]
    muts = l[5]
    linkerLen = len(linker)
    mLinkerLen = l[3] - l[2]
    mReadLen = rEnd - rStart
    if matches + muts >= linkerLen:  # mutations, deletions
        return l


## a mutation at the end is indistinguishable from a deletion we simply
## count it as a mutation because thats more likley
def correctLinkerLength(linker, rs, re, qs, qe, matches, mm):
    if mm == 1 and len(linker) - 1 == re - rs:
        return (rs, re + 1, qs, qe, matches, mm)
    else:
        return (rs, re, qs, qe, matches, mm)

def positionLinker(readPart, ll1, ll2, llinker):
    full = findLinker(readPart, llinker, 0.05)
    if full is not None:
        l1 = findLinker(readPart, ll1, 0.1)
        l2 = findLinker(readPart, ll2, 0.1)
        l1 = correctLinkerLength(ll1, *l1)
        l2 = correctLinkerLength(ll2, *l2)
        if l1 is not None and l2 is not None:
            if l1[1] + 6 == l2[0]:
                return (l1[0], l1[1], l2[0], l2[1])
            else:
                EXTRACTCOUNTER["L1not6fromL2"] += 1
        else:
            EXTRACTCOUNTER["noL1orL2"] += 1
    else:
        EXTRACTCOUNTER["noLL"] += 1


def getUMI(read, linker2End):
    acg = read[(linker2End + 6):(linker2End + 9)]
    gac = read[(linker2End + 17):(linker2End + 20)]
    if correctPrime5(acg, gac):
        umi = read[(linker2End + 9):(linker2End + 17)]
        return umi
    else:
        EXTRACTCOUNTER["noACGorGAC"] += 1


def getBCs(read, linkers):
    linker1Start = linkers[0] + PS
    linker1End = linkers[1] + PS
    linker2Start = linkers[2] + PS
    linker2End = linkers[3] + PS
    umi = getUMI(read, linker2End)
    if umi is not None:
        bc1 = read[(linker1Start - 6):linker1Start]
        bc2 = read[linker1End:linker2Start]
        bc3 = read[linker2End:(linker2End + 6)]
        return (bc1, bc2, bc3, umi)


def correctBC(bc):
    if BARCODESMAP.get(bc) is not None:
        BCCOUNTER['exact'] += 1
        return bc
    distances = filter(lambda m: distance.levenshtein(bc, m) <= 1, BARCODES)
    if len(distances) == 1:
        BCCOUNTER['mm1'] += 1
        return distances[0]


def correctPrime5(pre, suff):
    dist = distance.levenshtein(PRIME5, pre+suff)
    PRIME5COUNTER[dist] += 1
    if dist <= 2:
        return True
    else:
        return False

def extractBCsUMI(read, ll1, ll2, llinker):
    readPart = read[PS:PE]
    linkers = positionLinker(readPart, ll1, ll2, llinker)
    if linkers is not None and len(linkers) == 4:
        bcsU = getBCs(read, linkers)
        if bcsU is not None:
           bc1, bc2, bc3, umi = correctBC(bcsU[0]), correctBC(bcsU[1]), correctBC(bcsU[2]), bcsU[3]
           EXTRACTCOUNTER["extracted"] += 1
           return formatBC(bc1, bc2, bc3, umi)

def formatBC(bc1, bc2, bc3, umi):
    if bc1 is not None and bc2 is not None and bc3 is not None and umi is not None:
       return "BC:"+bc1+bc2+bc3+":UMI:"+umi
    
def writeCounter(counter):
    result = ""
    for k,v in sorted(counter.items()):
        result += str(k)+"\t"+str(v)+"\n"
    return result





