import os
import distance

from cutadapt.align import Aligner
from cutadapt.adapters import ANYWHERE


PS = 5
PE = 48
L1A = "TAGCCATCGCATTGC"
L2A = "TACCTCTGAGCTGAA"
L2B = "TACCTCTTATCTCTT"
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
STATSCOUNTER = {"noLL": 0, "noL1orL2": 0,
                "L1not6fromL2": 0, "noACGorGAC": 0, "extracted": 0}

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

# returns None or (ls,le) => for stats
# tuple: (l1s,l1e,l2s,l2e)
# check l1e + 6 == l2s


def positionLinker(readPart, ll1, ll2, llinker):
    full = findLinker(readPart, llinker, 0.05)
    if full is not None:
        l1 = findLinker(readPart, ll1, 0.1)
        l2 = findLinker(readPart, ll2, 0.1)
        if l1 is not None and l2 is not None:
            if l1[1] + 6 == l2[0]:
                return (l1[0], l1[1], l2[0], l2[1])
            else:
                STATSCOUNTER["L1not6fromL2"] += 1
        else:
            STATSCOUNTER["noL1orL2"] += 1
    else:
        STATSCOUNTER["noLL"] += 1


def getUMI(read, linker2End):
    acg = read[(linker2End + 6):(linker2End + 9)]
    gac = read[(linker2End + 17):(linker2End + 20)]
    if acg == "ACG" and gac == "GAC":
        umi = read[(linker2End + 9):(linker2End + 17)]
        return umi
    else:
        STATSCOUNTER["noACGorGAC"] += 1


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
        return bc
    distances = filter(lambda m: distance.levenshtein(bc, m) <= 1, BARCODES)
    if len(distances) == 1:
        return distances[0]


def extractBCsUMI(read, ll1, ll2, llinker):
    readPart = read[PS:PE]
    linkers = positionLinker(readPart, ll1, ll2, llinker)
    if linkers is not None and len(linkers) == 4:
        bcsU = getBCs(read, linkers)
        if bcsU is not None:
           bc1, bc2, bc3, umi = correctBC(bcsU[0]), correctBC(bcsU[1]), correctBC(bcsU[2]), bcsU[3]
           STATSCOUNTER["extracted"] += 1
           return formatBC(bc1, bc2, bc3, umi)

def formatBC(bc1, bc2, bc3, umi):
    if bc1 is not None and bc2 is not None and bc3 is not None and umi is not None:
       return "BC:"+bc1+bc2+bc3+":UMI:"+umi


def writeStatsCounter():
    result = ""
    for k,v in sorted(STATSCOUNTER.items()):
        result += k+"\t"+str(v)+"\n"
    return result








