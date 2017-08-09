import pytest


from surecellumi.mergefq import *
import surecellumi.umi

def test_merge():
    ll1 = umi.getLinker1("L1A")
    ll2 = umi.getLinker2("L2A")
    llinker = umi.getLLinker("L1A", "L2A")
    mergeUmiFiles("tests/data/r1.fq", "tests/data/r2.fq", "tests/data/merged.fq", "tests/data/stats.txt", ll1, ll2, llinker)
    assert open("tests/data/merged.fq",'r').read() == open("tests/data/merged.fq.result",'r').read()
    assert open("tests/data/stats.txt",'r').read() == open("tests/data/stats.txt.result",'r').read()
