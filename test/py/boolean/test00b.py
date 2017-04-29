""" A Gentle introduction to first steps of 3D Merge algorithm. Part 2. """
from larlib import *

V,[VV,EV,FV,CV] = larCuboids([2,2,1],True)
cubeGrid = Struct([(V,FV,EV)],"cubeGrid")
#cubeGrids = Struct(2*[cubeGrid,t(.5,.0,.5),r(0,0,PI/6)])
cubeGrids = Struct(2*[cubeGrid,t(.5,.5,.5),r(0,0,PI/6)])

V,FV,EV = struct2lar(cubeGrids)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV))))
W,FW,EW = partition(V,FV,EV)

WW = AA(LIST)(range(len(W)))
submodel = STRUCT(MKPOLS((W,EW)))
VIEW(larModelNumbering(1,1,1)(W,[WW,EW,FW],submodel,0.6)) 

SB_2 = SBoundary2(EW,FW)
SB_2.todense()

FE = [list(SB_2[:,f].tocoo().row) for f in range(SB_2.shape[1])]
triangleSet = boundaryTriangulation(W,FW,EW,FE)
TW,FT = triangleIndices(triangleSet,W)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,TW))))

SB_3 = SBoundary3(W,EW,FW)
print SB_3.todense()

