import pytest


from surecellumi.mergefq import *


def test_merge():
    mergeUmiFiles("tests/data/r1.fq", "tests/data/r2.fq", "tests/data/merged.fq", "tests/data/stats.txt")
    assert open("tests/data/merged.fq",'r').read() == open("tests/data/merged.fq.result",'r').read()
    assert open("tests/data/stats.txt",'r').read() == open("tests/data/stats.txt.result",'r').read()
