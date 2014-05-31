""" Computation of the boundary of a simplicial grid """
import sys; sys.path.insert(0, 'lib/py/')
from larcc import *
from largrid import *

V,CV = larSimplexGrid1((2,2,2))
bases = larSimplicialStack(CV)
VV,EV,FV,CV = bases
boundaryCells = signedBoundaryCells(V,CV,FV)

def swap(mylist): return [mylist[1]]+[mylist[0]]+mylist[2:]
orientedBoundary = [FV[-k] if k<0 else swap(FV[k]) for k in boundaryCells]
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,orientedBoundary))))
submodel = EXPLODE(1.05,1.05,1.05)(MKPOLS((V,orientedBoundary)))
VIEW(larModelNumbering(V,bases,submodel))
