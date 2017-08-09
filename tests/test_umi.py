import pytest

import distance
from cutadapt.align import Aligner
from cutadapt.adapters import ANYWHERE

from surecellumi.umi import * 

l1 = getLinker1("L1A")
l2 = getLinker2("L2A")
llinker = getLLinker("L1A", "L2A")

def test_findLinker():
    r1 = "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"
    #"GGGGGG"+"TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"+"GGGGGG"
    perfect = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    llinker = getLLinker("L1A", "L2A")
    assert findLinker(perfect, llinker, 0.05) == (6, 42, 0, 36, 36, 0)
    mut1 = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGTTGAA" + "GGGGGG"  # C ==> T
    assert findLinker(mut1, llinker, 0.05) == (6, 42, 0, 36, 35, 1)
    mut2 = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGTAGAA" + "GGGGGG"  # CT ==> TA
    assert findLinker(mut2, llinker, 0.05) == None
    mutStart = "GGGGGG" + "CAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(mutStart, llinker, 0.05) == (6, 42, 0, 36, 35, 1)
    mutEnd = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAT" + "GGGGGG"  # might have to search for next ATG to isolate BC3
    assert findLinker(mutEnd, llinker, 0.05) == (6, 41, 0, 36, 35, 1)
    delStart = "GGGGGG" + "AGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(delStart, llinker, 0.05) == (5, 41, 0, 36, 35, 1)
    delEnd = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGA" + "GGGGGG"
    assert findLinker(delEnd, llinker, 0.05) == (6, 41, 0, 36, 35, 1)
    insert = "GGGGGG" + "TAGGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"  # one G at start
    assert findLinker(insert, llinker, 0.05) == (6, 43, 0, 36, 36, 1)
    perfect = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(perfect, llinker, 0.05) == (6, 42, 0, 36, 36, 0)
    mutEnd1 = "GGGGGG" + "TAGCCATCGCATTGAAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(mutEnd1, llinker, 0.05) == (6, 42, 0, 36, 35, 1)
    perfect = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(perfect, llinker, 0.05) == (6, 42, 0, 36, 36, 0)


def test_positionLinker():
        # NNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTTTT
    r1 = "AAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT"
    assert positionLinker(r1[PS:PE], l1, l2, llinker) == (6, 21, 27, 42)
    r2 = "AAAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT"
    assert positionLinker(r2[PS:PE], l1, l2, llinker) == (7, 22, 28, 43)
    r3 = "AAAAAAATAGCCATCGCATTGGGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT" #mutation at end of linker1
    assert positionLinker(r3[PS:PE], l1, l2, llinker) == (7, 22, 28, 43)
    r4 = "AAAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGATCCCCCCACGCCCCCCCCGACTTTTTTTTTTTT" #mutation at end of linker2
    assert positionLinker(r4[PS:PE], l1, l2, llinker) == (7, 22, 28, 43)
    
# ANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTTT
# CTNNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTT
# GCANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTT
# TGCGNNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTT
# ATCGANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTT


def test_getBCs():
    bcm = ("AAAAAA", "GGGGGG", "CCCCCC", "CCCCCCCC")
    r1 = "AAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT"
    assert getBCs(r1, (6, 21, 27, 42)) == bcm
    assert getBCs("GG" + r1, (8, 23, 29, 44)) == bcm
   
def test_correctBC():
    assert "CGGTCC" == correctBC("CGGTCC")
    assert "CGGTCC" == correctBC("CGGACC")
    assert None == correctBC("CGAACC")

def test_data():
                                #     TACCTCTTATCTCTT???
         #ANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTTT
#    r1 = "AATCCGGTAGCCATCGCATTGCTTCTTGTACCTCTTATCTCTTACCCACACCTCCCACCTCACTTTTT"
#    assert positionLinker(r1[PS:PE], l1, getLinker2("L2B"), getLLinker("L1A","L2B")) == (2, 17, 23, 38)    
#    assert getBCs(r1, (2, 17, 23, 38)) == ('ATCCGG', 'TTCTTG', 'ACCCAC', 'TCCCACCT')
             #  AAGTAGCCATCGCATTGCAAGCCATACCTCTGAGCTGAAGCTC 
          #ANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTTT
#    r2 = "CTCGAAAGTAGCCATCGCATTGCAAGCCATACCTCTGAGCTGAAGCTCCCACGTCACCACAGACTTTT"
#    assert positionLinker(r2[PS:PE], l1, getLinker2("L2A"), getLLinker("L1A","L2A")) == (3, 18, 24, 39)    
#    assert getBCs(r2, (3, 18, 24, 39)) == ('CGAAAG', 'AAGCCA', 'GCTCCC', 'TCACCACA')
#    assert extractBCsUMI(r2, l1, getLinker2("L2A"), getLLinker("L1A","L2A")) == "BC:CGAAAGAAGCCAGCTCCC:UMI:TCACCACA"
    r3 = "AGGTTAGCCATCGCATTGCCTTACGTTACCTCTGAGCTGAAGGATTGACGCTGCTCATGACTTTTTTT" # too far left, L1 < 7nt > L2
    assert positionLinker(r3[PS:PE], l1, getLinker2("L2A"), getLLinker("L1A","L2A")) == None


