import pytest

import distance
from cutadapt.align import Aligner
from cutadapt.adapters import ANYWHERE

from surecellumi.umi import * 

def test_findLinker():
    r1 = "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"
    #"GGGGGG"+"TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA"+"GGGGGG"
    perfect = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(perfect, LINKER, 0.05) == (6, 42, 0, 36, 36, 0)
    mut1 = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGTTGAA" + "GGGGGG"  # C ==> T
    assert findLinker(mut1, LINKER, 0.05) == (6, 42, 0, 36, 35, 1)
    mut2 = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGTAGAA" + "GGGGGG"  # CT ==> TA
    assert findLinker(mut2, LINKER, 0.05) == None
    mutStart = "GGGGGG" + "CAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(mutStart, LINKER, 0.05) == (6, 42, 0, 36, 35, 1)
    mutEnd = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAT" + "GGGGGG"  # might have to search for next ATG to isolate BC3
    assert findLinker(mutEnd, LINKER, 0.05) == (6, 41, 0, 36, 35, 1)
    delStart = "GGGGGG" + "AGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(delStart, LINKER, 0.05) == (5, 41, 0, 36, 35, 1)
    delEnd = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGA" + "GGGGGG"
    assert findLinker(delEnd, LINKER, 0.05) == (6, 41, 0, 36, 35, 1)
    insert = "GGGGGG" + "TAGGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"  # one G at start
    assert findLinker(insert, LINKER, 0.05) == (6, 43, 0, 36, 36, 1)
    perfect = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(perfect, LINKER, 0.05) == (6, 42, 0, 36, 36, 0)
    perfect = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(perfect, LINKER, 0.05) == (6, 42, 0, 36, 36, 0)
    perfect = "GGGGGG" + "TAGCCATCGCATTGCAGCAGCTACCTCTGAGCTGAA" + "GGGGGG"
    assert findLinker(perfect, LINKER, 0.05) == (6, 42, 0, 36, 36, 0)


def test_positionLinker():
        # NNNNNNTAGCCATCGCATTGCNNNNNNTACCTCTGAGCTGAANNNNNNACGNNNNNNNNGACTTTTTTTTTTTT
    r1 = "AAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT"
    assert positionLinker(r1[PS:PE]) == (1, 16, 22, 37)
    r2 = "AAAAAAATAGCCATCGCATTGCGGGGGGTACCTCTGAGCTGAACCCCCCACGCCCCCCCCGACTTTTTTTTTTTT"
    assert positionLinker(r2[PS:PE]) == (2, 17, 23, 38)
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
