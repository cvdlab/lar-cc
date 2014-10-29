""" Oriented cuboidal cells """
""" Oriented cuboidal cells """
import sys;sys.path.insert(0, 'lib/py/')
from larcc import *



def orientedBoundaryCells(V,(VV,EV,FV,CV)):
    boundaryMat = signedCellularBoundary(V,[VV,EV,FV,CV])
    chainCoords = csc_matrix((len(CV), 1))
    for cell in range(len(CV)): chainCoords[cell,0] = 1
    boundaryCells = list((boundaryMat * chainCoords).tocoo().row)
    orientations = list((boundaryMat * chainCoords).tocoo().data)
    return zip(orientations,boundaryCells)

def normalVector(V,facet):
    v0,v1,v2 = facet[:3]
    return VECTPROD([ DIFF([V[v1],V[v0]]), DIFF([V[v2],V[v0]]) ])



# cuboidal grid
V,bases = larCuboids([5,5,3],True)
[VV,EV,FV,CV] = bases
boundaryCellspairs = orientedBoundaryCells(V,[VV,EV,FV,CV])

orientedBoundary = [FV[face] if sign>0 else swap(FV[face]) for (sign,face) in boundaryCellspairs]
normals = [ normalVector(V,facet)  for facet in orientedBoundary ]
facetCentroids = [CCOMB([V[v] for v in facet]) for facet in orientedBoundary]
appliedNormals = [[centroid,SUM([centroid,normal])] for (centroid,normal) in zip(facetCentroids,normals)]
normalVectors = AA(POLYLINE)(appliedNormals)

orientedQuads = [[sign,FV[face]] if sign>0 else [sign,swap(FV[face])] for (sign,face) in boundaryCellspairs]
FVtriangles = CAT([[[v0,v1,v2],[v2,v1,v3]] if sign>0 else [[v0,v1,v2],[v0,v2,v3]]
            for (sign,[v0,v1,v2,v3]) in orientedQuads])

VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FVtriangles))+normalVectors))

VIEW(MKPOL([V, [[v+1 for v in face] for face in FVtriangles], None ]))
