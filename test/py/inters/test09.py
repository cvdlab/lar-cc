""" Biconnected components from orthogonal LAR model """
import sys
sys.path.insert(0, 'lib/py/')
from inters import *
colors = [CYAN, MAGENTA, YELLOW, RED, GREEN, ORANGE, PURPLE, WHITE, BLACK, BLUE]

lines = randomLines(800,0.2)
V,EV = lines2lar(lines)
model = V,EV

VV = vertices2vertices(model)
leaves = [k for k,vv in enumerate(VV) if len(vv)==1]
EV_ = [[v1,v2]  for v1,v2 in EV if set(leaves).intersection([v1,v2]) == set()]
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,EV_))))
EV = AA(sorted)(EV_)

from bool1 import larRemoveVertices
V,EV = larRemoveVertices(V,EV)
VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV],submodel,0.015))
