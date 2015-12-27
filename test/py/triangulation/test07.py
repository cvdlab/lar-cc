
from larlib import *

""" random 1-boundary generation """
from larlib import *

import sys
sys.path.insert(0, '/Users/paoluzzi/Documents/dev/lar-cc/test/py/larcc/')
from test16 import *

EV = AA(list)(cells)
V,FV,EV = larPair2Triple((V,EV))

bcycles,bverts = boundaryCycles(range(len(EV)),EV)
VIEW(STRUCT(AA(POLYLINE)([[V[v] for v in verts+[verts[0]]] for verts in bverts])))

colors = [CYAN,MAGENTA,WHITE,RED,YELLOW,GRAY,GREEN,ORANGE,BLUE,PURPLE,BROWN,BLACK]
components = [COLOR(colors[k%12])(face) for k,face in enumerate(MKFACES((V,FV,EV)))]
VIEW(STRUCT(components))

