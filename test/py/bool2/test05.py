""" non-valid -> valid solid representation of a space partition """

""" Two unit cubes """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *
glass = MATERIAL([1,0,0,0.1,  0,1,0,0.1,  0,0,1,0.1, 0,0,0,0.1, 100])

V,[VV,EV,FV,CV] = larCuboids([1,1,1],True)
cube1 = Struct([(V,FV,EV)],"cube1")
#twoCubes = Struct([cube1,t(-1,.5,1),cube1])     # other test example
twoCubes = Struct([cube1,t(.5,.5,.5),cube1])
#twoCubes = Struct([cube1,t(.5,.5,0),cube1])    # other test example
#twoCubes = Struct([cube1,t(.5,0,0),cube1])        # other test example
V,FV,EV = struct2lar(twoCubes)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV))))

quadArray = [[V[v] for v in face] for face in FV]
boxes = containmentBoxes(quadArray)
hexas = AA(box2exa)(boxes)
parts = boxBuckets(boxes)

    
W,FW,EW = spacePartition(V,FV,EV, parts)

from architectural import *
polylines = lar2polylines((W,FW))
VIEW(EXPLODE(1.2,1.2,1.2)(AA(POLYLINE)(polylines)))

WW = AA(LIST)(range(len(W)))
submodel = STRUCT(MKPOLS((W,EW)))
VIEW(larModelNumbering(1,1,1)(W,[WW,EW,FW],submodel,0.5))
