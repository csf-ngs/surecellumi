#!/usr/bin/env python

import gzip, logging
import umi
from itertools import izip



def op(path):
    if path.endswith("gz"):
       return gzip.open(path, "r")
    else:
       return open(path, "r")

def mergeUmiFiles(read1file, read2file, outmerge, outstats, outbc, outprime5, ll1, ll2, llinker):
    logging.info("extracting umis with: "+ll1+" "+ll2+" "+llinker)
    lineNr = 0
    with op(read1file) as f1, op(read2file) as f2, open(outmerge, "w") as o:
        ignore = f1.next()
        for l1,l2 in izip(f1,f2):
           if lineNr % 1e6 == 0:
              logging.info("extracting umi lineNr: "+str(lineNr))
           if lineNr % 4 == 0:
              try:
                 bcu = umi.extractBCsUMI(l1.strip(), ll1, ll2, llinker)
              except:
                 print("error in line: "+str(lineNr)+" "+read1file+" extract UMI\n"+l1.strip())
                 raise
              if bcu is not None:
                  l2 = l2.strip() + ":" + bcu + "\n"
              else:
                  l2 = l2.strip() + ":" + l1
           lineNr += 1
           o.write(l2)
        o.write(f2.next())
    o.close()
    logging.info("extracted last umi lineNr: "+str(lineNr+1))
    logging.info("extract stats:\n"+umi.writeCounter(umi.EXTRACTCOUNTER))
    logging.info("bc stats:\n"+umi.writeCounter(umi.BCCOUNTER))
    logging.info("prime5 stats:\n"+umi.writeCounter(umi.PRIME5COUNTER))
    with open(outstats, "w") as o:
         o.write(umi.writeCounter(umi.EXTRACTCOUNTER))
    with open(outbc, "w") as o:
         o.write(umi.writeCounter(umi.BCCOUNTER))
    with open(outprime5, "w") as o:
         o.write(umi.writeCounter(umi.PRIME5COUNTER))
    o.close()

