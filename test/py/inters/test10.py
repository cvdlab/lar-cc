""" Biconnected components from orthogonal LAR model """
import sys
sys.path.insert(0, 'lib/py/')
from inters import *
from bool1 import larRemoveVertices
from hospital import surfIntegration
from iot3d import polyline2lar

filename = "test/py/inters/test1.svg"
lines = svg2lines(filename)
VIEW(STRUCT(AA(POLYLINE)(lines)))

V,EV = lines2lar(lines)
print "\nV =",V
print "\nEV =",EV
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,EV))))
model = V,EV
VV = vertices2vertices(model)
leaves = [k for k,vv in enumerate(VV) if len(vv)==1]
EV_ = [[v1,v2]  for v1,v2 in EV if set(leaves).intersection([v1,v2]) == set()]
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,EV_))))

EV = list(set(AA(tuple)(sorted(AA(sorted)(EV_))))) 
V,EV = larRemoveVertices(V,EV)
VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV],submodel,0.10))

model = V,EV
FV = facesFromComponents((V,EV))
areas = surfIntegration((V,FV,EV))
boundaryArea = max(areas)
faces = [FV[f] for f,area in enumerate(areas) if area!=boundaryArea]
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,faces+EV)) + AA(MK)(V)))

V,FV,EV = polyline2lar([[ V[v] for v in FV[areas.index(boundaryArea)] ]])
VIEW(STRUCT(MKPOLS((V,EV))))
