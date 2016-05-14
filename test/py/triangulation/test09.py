""" Test example of LAR of a 2-complex with non-contractible and non-manifold cells"""
from larlib import *

filename = "test/svg/inters/test0.svg"
lines = inters.svg2lines(filename)
V,FV,EV,polygons = larFromLines(lines,True)

VV = AA(LIST)(range(len(V)))
hpc = STRUCT(MKPOLS((V,EV)))        
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],hpc,0.2)) 
