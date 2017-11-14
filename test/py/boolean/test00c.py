""" A Gentle introduction to first steps of 3D Merge algorithm. Part 2. """
from larlib import * 

V,[VV,EV,FV,CV] = larCuboids([1,1,1],True)
cubeGrid = Struct([(V,FV,EV)],"cubeGrid")
cubeGrids = Struct(2*[cubeGrid,t(.5,.5,.5),r(0,0,PI/6)])#,r(0,PI/3,0)])
#cubeGrids = Struct(2*[cubeGrid,t(.5,.5,.0)])

V,FV,EV = struct2lar(cubeGrids)
VV = AA(LIST)(range(len(V)))
VIEW(STRUCT(MKPOLS((V,[range(8),range(8,16)]))))
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV+EV+VV))))
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

B_3 = Boundary3(W,EW,FW)
print B_3.todense()

hpcCells = MKSOLID(W,FW,EW)
VIEW(EXPLODE(1.5,1.5,1.5)(hpcCells))

FC = B_3.todense()
cell2faces = TRANS(FC.tolist())
CV = [set(CAT([FW[f] for (f,val) in enumerate(cell) if val!=0])) for cell in cell2faces]

