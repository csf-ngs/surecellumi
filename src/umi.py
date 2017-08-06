import unittest

import distance
from cutadapt.align import Aligner
from cutadapt.adapters import ANYWHERE

PS=5
PE=48
L1="TAGCCATCGCATTGC" 
L2="TACCTCTGAGCTGAA" 
LINKER = L1+"NNNNNN"+L2 

def loadBarcodes():
    barcodes = []
    with open("barcodes.txt") as b:
         for l in b:
             barcodes.append(l.strip())    
    return barcodes 

BARCODES = loadBarcodes()
BARCODESMAP = dict(zip(BARCODES,BARCODES))
STATSCOUNTER = { "noLL": 0, "noL1orL2": 0, "L1not6fromL2": 0, "noACGorGAC": 0, "extracted": 0 }

##mutation rate full linker: length is 6 + 15+6+15 + 6 = 48: 0.05
##mutation rate linker part: 0.1 
def findLinker(readPart, linker, mut):
    aligner = Aligner(readPart, mut, ANYWHERE, wildcard_ref=True,wildcard_query=True)
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
    if matches + muts >= linkerLen: #mutations, deletions
       return l
  
##returns None or (ls,le) => for stats
##tuple: (l1s,l1e,l2s,l2e) 
##check l1e + 6 == l2s
def positionLinker(readPart):
    full = findLinker(readPart, LINKER, 0.05)
    if full is not None:
       l1 = findLinker(readPart, L1, 0.1)
       l2 = findLinker(readPart, L2, 0.1)
       if l1 is not None and l2 is not None:
          if l1[1] + 6 == l2[0]:
             return (l1[0],l1[1],l2[0],l2[1])
          else:
             STATSCOUNTER["L1not6fromL2"] += 1
       else:
          STATSCOUNTER["noL1orL2"] += 1
    else:
       STATSCOUNTER["noLL"] += 1

def getUMI(read, linker2End):
    acg = read[(linker2End+6):(linker2End+9)]
    gac = read[(linker2End+17):(linker2End+20)]
    if acg == "ACG" and gac == "GAC":
       umi = read[(linker2End+9):(linker2End+17)]    
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
       bc1 = read[(linker1Start-6):linker1Start]
       bc2 = read[linker1End:linker2Start]
       bc3 = read[linker2End:(linker2End+6)]
       return (bc1, bc2, bc3, umi)

def correctBC(bc):
    if BARCODESMAP.get(bc) is not None:
       return bc
    distances = filter(lambda m: distance.levenshtein(bc, m) <=  1, BARCODES)      
    if len(distances) == 1:
       return distances[0]

def extractBCsUMI(read):
    readPart = read[PS:PE]
    linkers = positionLinker(readPart)
    if linkers is not None and len(linkers) == 4:
       bcsU = getBCs(read, linkers) 
       bc1, bc2, bc3, umi = correctBC(bcsU[0]), correctBC(bcsU[1]), correctBC(bcsU[2]), bcsU[3] 
       STATSCOUNTER["extracted"] += 1
       return (bc1, bc2, bc3, umi)        

class TestLinkerFinder(unittest.TestCase):
     
     def test_findLinker(self):
	 r1 = "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"
	 #"GGGGGG"+"TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"+"GGGGGG"
	 perfect =  "GGGGGG"+"TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"+"GGGGGG"
	 self.assertEqual(findLinker(perfect, LINKER, 0.05), (6,42,0,36,36,0))
	 mut1 =     "GGGGGG"+"TAGCCATCGCATTGCAGCAGCTACCTCTGAGTTGAA"+"GGGGGG" #C => T
	 self.assertEqual(findLinker(mut1, LINKER, 0.05), (6,42,0,36,35,1))
         mut2 =     "GGGGGG"+"TAGCCATCGCATTGCAGCAGCTACCTCTGAGTAGAA"+"GGGGGG" #CT => TA
         self.assertEqual(findLinker(mut2, LINKER, 0.05), None)
         mutStart =  "GGGGGG"+"CAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"+"GGGGGG" 
         self.assertEqual(findLinker(mutStart, LINKER, 0.05), (6,42,0,36,35,1))
         mutEnd =  "GGGGGG"+"TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAT"+"GGGGGG" #might have to search for next ATG to isolate BC3
         self.assertEqual(findLinker(mutEnd, LINKER, 0.05), (6,41,0,36,35,1))
         delStart =  "GGGGGG"+"AGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"+"GGGGGG" 
         self.assertEqual(findLinker(delStart, LINKER, 0.05), (5,41,0,36,35,1))
         delEnd =  "GGGGGG"+"TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGA"+"GGGGGG"
         self.assertEqual(findLinker(delEnd, LINKER, 0.05), (6,41,0,36,35,1))
         insert =  "GGGGGG"+"TAGGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"+"GGGGGG" #one G at start
         self.assertEqual(findLinker(insert, LINKER, 0.05), (6,43,0,36,36,1)) 
         perfect =  "GGGGGG"+"TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"+"GGGGGG"
         self.assertEqual(findLinker(perfect, LINKER, 0.05), (6,42,0,36,36,0))
         perfect =  "GGGGGG"+"TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"+"GGGGGG"
         self.assertEqual(findLinker(perfect, LINKER, 0.05), (6,42,0,36,36,0))
         perfect =  "GGGGGG"+"TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"+"GGGGGG"
         self.assertEqual(findLinker(perfect, LINKER, 0.05), (6,42,0,36,36,0))

     def test_positionLinker(self):
              #NNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTTTT
         r1 = "AAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT"
         self.assertEqual(positionLinker(r1[PS:PE]), (1, 16, 22, 37))
         r2 = "AAAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT"  
         self.assertEqual(positionLinker(r2[PS:PE]), (2, 17, 23, 38))
#ANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTTT
#CTNNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTT
#GCANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTT
#TGCGNNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTT
#ATCGANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTT     
 
     def test_getBCs(self):
         bcm = ("AAAAAA","GGGGGG","CCCCCC","CCCCCCCC")
         r1 = "AAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT"
         self.assertEqual( getBCs(r1, (1, 16, 22, 37)), bcm )
         self.assertEqual( getBCs("GG"+r1, (3, 18, 24, 39)), bcm )

     def test_correctBC(self):
         self.assertEqual( "CGGTCC", correctBC("CGGTCC"))
         self.assertEqual( "CGGTCC", correctBC("CGGACC"))
         self.assertEqual( None    , correctBC("CGAACC"))

     def nohting(self):	 
	 mut1 = "GGGGGG"++"GGGGGG"
         self.assertEqual(findLinker(mut1, l1, 0.1), (3,18)) 
	 mutStart = "GGGGGG"+"CAGCCATCGCATTGC"+"GGGGGG"
	 self.assertEqual(findLinker(mutStart, l1, 0.1), (3,18))
         mutEnd1 = "GGGGGG"+"TAGCCATCGCATTGT"+"GGGGGG"
 	 self.assertEqual(findLinker(mutEnd1, l1, 0.1), (3,18))
         mutEnd2 = "GGGGGG"+"TAGCCATCGCATTGG"+"GGGGGG"      
         self.assertEqual(findLinker(mutEnd2, l1, 0.1), (3,18))
         insert =  "GGGGGG"+"TAGCCATCGGCATTGC"+"GGGGGG" 
         self.assertEqual(findLinker(insert, l1, 0.1), (3,18))
         deletion = "GGGGGG"+"TAGCCATGCATTGC"+"GGGGGG"
         self.assertEqual(findLinker(deletion, l1, 0.1), (3,18))
	 delStart = "GGGGGG"+"AGCCATCGCATTGC"+"GGGGGG"
         self.assertEqual(findLinker(delStart, l1, 0.1), (3,18))
         delEnd = "GGGGGG"+"TAGCCATCGCATTG"+"GGGGGG"
	 self.assertEqual(findLinker(delEnd, l1, 0.1), (3,18))
         mut2 = "GGGGGG"+"TAGGCATCGCAATGC"+"GGGGGG"
         self.assertEqual(findLinker(mut2, l1, 0.1), None)


if __name__ == '__main__':
     unittest.main()
