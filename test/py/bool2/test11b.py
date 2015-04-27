

""" Visualization of indices of the boundary triangulation """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *

V,[VV,EV,FV,CV] = larCuboids([2,2,1],True)
cube1 = Struct([(V,FV,EV)],"cube1")
twoCubes = Struct(2*[cube1,t(.5,.5,0)])

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
V,CV,FV,EV,CF,CE,COE = facesFromComponents(theModel)

print ">> 10: CF",[(k,cf) for k,cf in enumerate(CF)]

CF = sorted(list(set(AA(tuple)(AA(sorted)(CF)))))
cellLengths = AA(len)(CF)
boundaryPosition = cellLengths.index(max(cellLengths))
BF = CF[boundaryPosition]
del CF[boundaryPosition]
VIEW(EXPLODE(1.2,1.2,1.2)( MKTRIANGLES(W,FW) ))
VIEW(EXPLODE(1.5,1.5,1.5)( MKTRIANGLES(V,[FV[f] for f in BF]) ))
VIEW(EXPLODE(2,2,2)([ MKSOLID(V,[FV[f] for f in cell]) for cell in CF]))

WW = AA(LIST)(range(len(W)))
submodel = SKEL_1(STRUCT(MKPOLS((W,EW))))
VIEW(larModelNumbering(1,1,1)(W,[WW,EW,FW],submodel,0.6))

