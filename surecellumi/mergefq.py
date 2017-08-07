#!/usr/bin/env python

import gzip
import umi
from itertools import izip

def op(path):
    if path.endswith("gz"):
       return gzip.open(path, "r")
    else:
       return open(path, "r")

def mergeUmiFiles(read1file, read2file, outmerge, outstats):
    with op(read1file) as f1, op(read2file) as f2, open(outmerge, "w") as o:
        lineNr = 0
        ignore = f1.next()
        for l1,l2 in izip(f1,f2):
           if lineNr % 4 == 0:
              bcu = umi.extract(l1.strip())
              if bcu is not None:
                  l2 = l2.strip() + ":" + bcu + "\n"
              else:
                  l2 = l2.strip() + ":" + l1
           lineNr += 1
           o.write(l2)
    o.close()
    with open(outstats, "w") as o:
         o.write(umi.writeStatsCounter())
    o.close()

