""" Test example of LAR of a 2-complex with non-contractible and non-manifold cells"""
from larlib import *

filename = "test/svg/inters/test0.svg"
lines = svg2lines(filename)
V,EV = lines2lar(lines)
V,FV,EV = larPair2Triple((V,EV))
VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))        
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,0.3)) 
