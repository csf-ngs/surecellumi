import pytest


from surecellumi.mergefq import *
import surecellumi.umi

def test_merge():
    ll1 = umi.getLinker1("L1A")
    ll2 = umi.getLinker2("L2A")
    llinker = umi.getLLinker("L1A", "L2A")
    mergeUmiFiles("tests/data/r1.fq", "tests/data/r2.fq", "tests/data/merged.fq", "tests/data/extractstat.txt", "tests/data/bcstat.txt", ll1, ll2, llinker)
    assert open("tests/data/merged.fq",'r').read() == open("tests/data/merged.fq.result",'r').read()
    assert open("tests/data/extractstat.txt",'r').read() == open("tests/data/extractstat.txt.result",'r').read()
    assert open("tests/data/bcstat.txt",'r').read() == open("tests/data/bcstat.txt.result",'r').read()
