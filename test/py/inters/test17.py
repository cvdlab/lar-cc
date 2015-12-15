""" Generating the LAR of a set of non-intersecting cycles """
from larlib import *

sys.path.insert(0, 'test/py/inters/')
from test16 import *

lar = (V,FV,EV)

bcycles,_ = boundaryCycles(range(len(EV)),EV)
polylines = [[V[EV[e][1]] if e>0 else V[EV[-e][0]] for e in cycle ] for cycle in bcycles]
polygons = [polyline + [polyline[0]] for polyline in polylines]

complex = SOLIDIFY(STRUCT(AA(POLYLINE)(polygons)))
csrBoundaryMat = boundary(FV,EV)
for chain in chains:
    chainBoundary = COLOR(RED)(STRUCT(MKPOLS((V,[EV[e] 
                        for e in chain2BoundaryChain(csrBoundaryMat)(chain)]))))
    VIEW(STRUCT([complex, chainBoundary]))
