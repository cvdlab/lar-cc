""" Testing signed 2-boundary """
from larlib import *

lines = svg2lines("test/svg2lines/test.svg")
V,FV,EV,polygons = larFromLines(lines,True)
VIEW(STRUCT(MKTRIANGLES((V,FV,EV),color=True)))
submodel = mkSignedEdges((V,EV))
VV = AA(LIST)(range(len(V)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,0.25))

B = larSignedBoundary2(V,FV,EV)
for k in range(B.shape[0]):
    print k,B.todense()[k]

VIEW(STRUCT(MKTRIANGLES((V,FV[1:3],EV),color=True)))
