
from larlib import *

""" Visualization of indices of the boundary triangulation """

V,[VV,EV,FV,CV] = larCuboids([2,2,1],True)
"""
BF = [FV[f] for f in boundaryCells(CV,FV)]
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,BF))))
_,BE = larFacets((V,BF),dim=2)
cubeGrid = Struct([(V,BF,BE)],"cubeGrid")
"""
cubeGrid = Struct([(V,FV,EV)],"cubeGrid")
cubeGrids = Struct(3*[cubeGrid,t(.5,.5,0),r(0,0,PI/6)])

V,FV,EV = struct2lar(cubeGrids)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV))))
V,CV,FV,EV,CF,CE,COE,FE = thePartition(V,FV,EV)

"""
boundaryCells_2 = signedBoundaryCells(V,CV,FV)
boundaryFV = [list(FV[-k]) if k<0 else swap(list(FV[k])) for k in boundaryCells_2]
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,boundaryFV))))
VIEW(EXPLODE(1.5,1.5,1.5)(MKTRIANGLES((V,boundaryFV,EV))))
"""

cellLengths = AA(len)(CF)
boundaryPosition = cellLengths.index(max(cellLengths))
BF = CF[boundaryPosition]; del CF[boundaryPosition]; del CE[boundaryPosition]
BE = list({e for f in BF for e in FE[f]})

#Volume((V,[FV[f] for f in CF[0]]))

VIEW(EXPLODE(1.2,1.2,1.2)( MKTRIANGLES((V,[FV[f] for f in BF],[EV[e] for e in BE])) ))
VIEW(EXPLODE(1.2,1.2,1.2)([ MKSOLID(V,[FV[f] for f in cell],[EV[e] for e in set(CAT([FE[f] for f in cell]))]) for cell in CF]))
VIEW(EXPLODE(1.2,1.2,1.2)([STRUCT(MKPOLS((V,[EV[e] for e in cell]))) for cell in CE]))

