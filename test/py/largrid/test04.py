""" Computation of the boundary of a simplicial grid """
import sys; sys.path.insert(0, 'lib/py/')
from larcc import *
from largrid import *

V,simplices = larSimplexGrid1((2,2))
bases = larSimplicialStack(simplices)
boundaryCells = signedBoundaryCells(V,bases[-1],bases[-2])

def swap(mylist): return [mylist[1]]+[mylist[0]]+mylist[2:]
orientedBoundary = [bases[-2][-k] if k<0 else swap(bases[-2][k]) 
                  for k in boundaryCells]
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,orientedBoundary))))
submodel = STRUCT(MKPOLS((V,orientedBoundary)))
VIEW(larModelNumbering(V,bases,submodel))
