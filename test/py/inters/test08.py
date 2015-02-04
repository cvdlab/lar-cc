""" Biconnected components from orthogonal LAR model """
import sys
sys.path.insert(0, 'lib/py/')
from inters import *
colors = [CYAN, MAGENTA, YELLOW, RED, GREEN, ORANGE, PURPLE, WHITE, BLACK, BLUE]

lines = randomLines(800,0.2)
V,EV = lines2lar(lines)
model = V,EV
EVs = biconnectedComponent(model)
HPCs = [STRUCT(MKPOLS((V,ev))) for ev in EVs if len(ev)>1]
verts = list(set(CAT(CAT([ev for ev in EVs if len(ev)>1]))))

sets = [COLOR(colors[k%10])(hpc) for k,hpc in enumerate(HPCs)]
VIEW(STRUCT(sets+AA(MK)([V[v] for v in verts])))
VIEW(STRUCT(AA(POLYLINE)(lines)))
