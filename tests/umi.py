import pytest

import distance
from cutadapt.align import Aligner
from cutadapt.adapters import ANYWHERE


    
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



