
from larlib import *

""" random 1-boundary generation """
from larlib import *

import sys
sys.path.insert(0, '/Users/paoluzzi/Documents/dev/lar-cc/test/py/larcc/')
from test16 import *

EV = AA(list)(cells)
V,EVs = biconnectedComponent((V,EV))
FV = AA(COMP([sorted,list,set,CAT]))(EVs)
FV = sorted( FV,key=len,reverse=True )
EVs = sorted( EVs,key=len,reverse=True )
W = [eval(vcode(v)) for v in V]
latticeArray = computeCycleLattice(W,EVs)

viewLarComplexChain((V,FV,EV))

"""
bcycles,bverts = boundaryCycles(range(len(EW)),EW)
VIEW(STRUCT(AA(POLYLINE)([[V[v] for v in verts] for verts in bverts])))
"""
    
colors = [CYAN, MAGENTA, WHITE, RED, YELLOW, GRAY, GREEN, ORANGE, BLACK, BLUE, 
         PURPLE, BROWN]
components = [COLOR(colors[k%12])(STRUCT(MKPOLS((V,ev)))) for k,ev in enumerate(EVs)]
VIEW(STRUCT(components))

