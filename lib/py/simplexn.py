# -*- coding: utf-8 -*-
"""Module for facet extraction, extrusion and simplicial grids"""
from pyplasm import *
from lar2psm import *
from scipy import *

VOID = V0,CV0 = [[]],[[0]]    # the empty simplicial model

def cumsum(iterable):
    # cumulative addition: list(cumsum(range(4))) => [0, 1, 3, 6]
    iterable = iter(iterable)
    s = iterable.next()
    yield s
    for c in iterable:
        s = s + c
        yield s

def larExtrude(model,pattern):
    V, FV = model
    d, m = len(FV[0]), len(pattern)
    coords = list(cumsum([0]+(AA(ABS)(pattern))))
    offset, outcells, rangelimit = len(V), [], d*m
    for cell in FV:
        tube = [v + k*offset for k in range(m+1) for v in cell]  
        cellTube = [tube[k:k+d+1] for k in range(rangelimit)]
        outcells += [reshape(cellTube, newshape=(m,d,d+1)).tolist()]      
        
    outcells = AA(CAT)(TRANS(outcells))
    cellGroups = [group for k,group in enumerate(outcells) if pattern[k]>0 ]
    outVertices = [v+[z] for z in coords for v in V]
    outModel = outVertices, CAT(cellGroups)
    return outModel

def larSimplexGrid(shape):
    model = VOID
    for item in shape:
        model = larExtrude(model,item*[1])
    return model

def larSimplexFacets(simplices):
    out = []
    d = len(simplices[0])
    for simplex in simplices:
        out += [simplex[0:k]+simplex[k+1:d] for k in range(d)]
    out = sorted(out)
    return [facet for k,facet in enumerate(out[:-1]) if out[k] != out[k+1]] \
      + [out[-1]] 

if __name__ == "__main__":
   # example 1
   V = [[0,0],[1,0],[2,0],[0,1],[1,1],[2,1],[0,2],[1,2],[2,2]]
   FV = [[0,1,3],[1,2,4],[2,4,5],[3,4,6],[4,6,7],[5,7,8]]
   model = larExtrude((V,FV),4*[1,2,-3])
   VIEW(EXPLODE(1,1,1.2)(MKPOLS(model)))
   
   # example 2
   model = larExtrude( VOID, 6*[1] )
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
   model = larExtrude( model, 6*[1] )
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
   model = larExtrude( model, 6*[1] )
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
   
   # example 3
   model = larExtrude( VOID, 10*[1,-1] )
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
   model = larExtrude( model, 10*[1] )
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
   
   grid_2d = larSimplexGrid([3,3])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(grid_2d)))
   
   grid_3d = larSimplexGrid([2,3,4])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(grid_3d)))
   
   V,CV = larSimplexGrid([1,1,1])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))
   SK2 = (V,larSimplexFacets(CV))
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK2)))
   SK1 = (V,larSimplexFacets(SK2[1]))
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK1)))
   
