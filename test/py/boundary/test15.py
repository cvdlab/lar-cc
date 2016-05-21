""" Example of Offset generation for the 2-faces of a 2D complex """
from larlib import *
lines = svg2lines("test/svg2lines/test.svg")
V,FV,EV,polygons = larFromLines(lines,True)
VIEW(STRUCT(MKTRIANGLES((V,FV,EV),color=True)))

submodel = mkSignedEdges((V,EV))
VV = AA(LIST)(range(len(V)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,0.25))

newEdges = larOffset2D((V,FV,EV),offset=0.01)
VIEW(STRUCT(MKPOLS((V,EV)) + AA(COLOR(YELLOW))(AA(POLYLINE)(newEdges))))
