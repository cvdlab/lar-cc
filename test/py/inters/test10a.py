""" Biconnected components from orthogonal LAR model """
from larlib import *

print "\n drag your SVG file to  http://cvdlab.github.io/svg2lines"
print "then save it to <path/filename>.lines"
filename = raw_input('filename =')

#filename = "test/svg/inters/plan.lines"
#filename = "test/py/inters/building.svg"
#filename = "test/py/inters/complex.svg"

lines = lines2lines(filename)
VIEW(STRUCT(AA(POLYLINE)(lines)))
    
V,FV,EV,polygons = larFromLines(lines)
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,FV[:-1]+EV)) + AA(MK)(V)))

VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV[:-1]],submodel,0.05))


verts,faces,edges = polyline2lar([[ V[v] for v in FV[-1] ]])
VIEW(STRUCT(MKPOLS((verts,edges))))
