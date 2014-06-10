""" Oriented cuboidal and simplicial cells (same algorithm) """
import sys;sys.path.insert(0, 'lib/py/')
from larcc import *

V,bases = larCuboids([5,5,3],True)
[VV,EV,FV,CV] = bases
orientedBoundary = signedCellularBoundaryCells(V,AA(AA(REVERSE))([VV,EV,FV,CV]))
VIEW(EXPLODE(1.25,1.25,1.25)(MKPOLS((V,orientedBoundary))))

V,CV = larSimplexGrid1([5,5,3])
FV = larSimplexFacets(CV)
EV = larSimplexFacets(FV)
VV = AA(LIST)(range(len(V)))
bases = [VV,EV,FV,CV]
orientedBoundary = signedCellularBoundaryCells(V,bases)
VIEW(EXPLODE(1.25,1.25,1.25)(MKPOLS((V,orientedBoundary))))
