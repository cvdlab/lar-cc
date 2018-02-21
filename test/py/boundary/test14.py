""" Testing the extraction of a boundary chain """
from larlib import *

lines = svg2lines("test/svg2lines/test.svg")
V,FV,EV,polygons = larFromLines(lines,True)
VIEW(STRUCT(MKTRIANGLES((V,FV,EV),color=True)))
submodel = mkSignedEdges((V,EV))
VV = AA(LIST)(range(len(V)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,0.25))

orientations,boundaryCells = larSignedBoundary2Cells(V,FV,EV)([1,2])
orientedBoundaryCells = [EV[e] if sign==1 else REVERSE(EV[e]) 
                                                for sign,e in zip(orientations,boundaryCells)]

VIEW(STRUCT(MKTRIANGLES((V,FV[1:3],EV),color=True)))
VIEW(mkSignedEdges((V,orientedBoundaryCells)))
