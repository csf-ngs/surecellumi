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
    perfect = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(perfect, llinker, 0.05) == (6, 42, 0, 36, 36, 0)
    perfect = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(perfect, llinker, 0.05) == (6, 42, 0, 36, 36, 0)


def test_positionLinker():
        # NNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTTTT
    r1 = "AAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT"
    
    assert positionLinker(r1[PS:PE], l1, l2, llinker) == (1, 16, 22, 37)
    r2 = "AAAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT"
    assert positionLinker(r2[PS:PE], l1, l2, llinker) == (2, 17, 23, 38)
# ANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTTT
# CTNNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTT
# GCANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTT
# TGCGNNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTT
# ATCGANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTT


def test_getBCs():
    bcm = ("AAAAAA", "GGGGGG", "CCCCCC", "CCCCCCCC")
    r1 = "AAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT"
    assert getBCs(r1, (1, 16, 22, 37)) == bcm
    assert getBCs("GG" + r1, (3, 18, 24, 39)) == bcm


def test_correctBC():
    assert "CGGTCC" == correctBC("CGGTCC")
    assert "CGGTCC" == correctBC("CGGACC")
    assert None == correctBC("CGAACC")

def test_data():
                                #     TACCTCTTATCTCTT???
         #ANNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTTT
    r1 = "AATCCGGTAGCCATCGCATTGCTTCTTGTACCTCTTATCTCTTACCCACACCTCCCACCTCACTTTTT"
    assert positionLinker(r1[PS:PE], l1, getLinker2("L2B"), getLLinker("L1A","L2B")) == (2, 17, 23, 38)    






