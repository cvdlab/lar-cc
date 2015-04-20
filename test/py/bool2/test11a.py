

""" Visualization of indices of the boundary triangulation """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *

V,[VV,EV,FV,CV] = larCuboids([1,1,1],True)
cube1 = Struct([(V,FV,EV)],"cube1")
twoCubes = Struct([cube1,t(.5,.5,.5),cube1])

glass = MATERIAL([1,0,0,0.1,  0,1,0,0.1,  0,0,1,0.1, 0,0,0,0.1, 100])

#twoCubes = Struct([cube1,t(-1,.5,1),cube1])     # other test example
#twoCubes = Struct([cube1,t(.5,.5,0),cube1])    # other test example
#twoCubes = Struct([cube1,t(.5,0,0),cube1])        # other test example

V,FV,EV = struct2lar(twoCubes)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV))))

quadArray = [[V[v] for v in face] for face in FV]
boxes = containmentBoxes(quadArray)
hexas = AA(box2exa)(boxes)
parts = boxBuckets(boxes)

W,FW,EW = spacePartition(V,FV,EV, parts)
WW = AA(LIST)(range(len(W)))
submodel = STRUCT(MKPOLS((W,EW)))
VIEW(larModelNumbering(1,1,1)(W,[WW,EW,FW],submodel,0.6))

VIEW(SKEL_1(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,FW)))))

theModel = W,FW,EW
V,CV,FV,EV,CF,CE = facesFromComponents(theModel)

VIEW(EXPLODE(1.2,1.2,1.2)( MKTRIANGLES(W,FW) ))
VIEW(EXPLODE(2,2,2)([ STRUCT(MKTRIANGLES(V,[FV[f] for f in cf])) for cf in [CF[0],CF[1],CF[3]]  ]))
VIEW(EXPLODE(1.2,1.2,1.2)( MKTRIANGLES(V,[FV[f] for f in CF[2]]) ))

WW = AA(LIST)(range(len(W)))
submodel = SKEL_1(STRUCT(MKPOLS((W,EW))))
VIEW(larModelNumbering(1,1,1)(W,[WW,EW,FW],submodel,0.6))

