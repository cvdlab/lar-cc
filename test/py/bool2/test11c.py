

""" Visualization of indices of the boundary triangulation """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *

V,[VV,EV,FV,CV] = larCuboids([1,1,1],True)
cube1 = Struct([(V,FV,EV)],"cube1")
twoCubes = Struct([cube1,t(.5,0,0),cube1])

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

triangleSets = boundaryTriangulation(V,FV)
VIEW(EXPLODE(1.2,1.2,1.2)([STRUCT([MKPOL([tria,[[1,2,3,4]],None]) for tria in triangleSet]) for triangleSet in triangleSets]))

CF = AA(list)(CF)
CE = AA(list)(CE)

VIEW(STRUCT(AA(STRUCT)(AA(MKPOLS)( DISTL([V,[[EV[c] for c in cell] for cell in CE[:-1] ]])))  ))
VIEW(EXPLODE(2,2,2) (AA(STRUCT)(AA(MKPOLS)( DISTL([V,[[EV[c] for c in cell] for cell in [CE[-1]] ]])))  ))
VIEW(EXPLODE(2,2,2) (AA(STRUCT)(AA(MKPOLS)( DISTL([V,[[FV[c] for c in cell] for cell in [CF[-1]] ]])))  ))

models = DISTL([V,[[FV[c] for c in cell] for cell in CF ]])
models = [boundaryTriangulation(*model) for model in models]

def MKCELL(model): 
    return STRUCT([ STRUCT([MKPOL([tria,[[1,2,3,4]],None]) for tria in triangleSet]) 
             for triangleSet in model ])

VIEW(EXPLODE(1.5,1.5,1.5)(AA(MKCELL)([model for model in models])))

WW = AA(LIST)(range(len(W)))
submodel = SKEL_1(STRUCT(MKPOLS((W,EW))))
VIEW(larModelNumbering(1,1,1)(W,[WW,EW,FW],submodel,0.6))

