""" Oriented cuboidal and simplicial cells (same algorithm) """
import sys;sys.path.insert(0, 'lib/py/')
from larcc import *

# cuboidal grid
V,bases = larCuboids([5,5,3],True)
[VV,EV,FV,CV] = bases
orientedBoundary = signedCellularBoundaryCells(V,AA(AA(REVERSE))([VV,EV,FV,CV]))
cells = [FV[f] if sign==1 else REVERSE(FV[f])  for (sign,f) in zip(*orientedBoundary)]
VIEW(EXPLODE(1.25,1.25,1.25)(MKPOLS((V,cells))))

# simplicial grid
V,CV = larSimplexGrid1([5,5,3])
FV = larSimplexFacets(CV)
EV = larSimplexFacets(FV)
VV = AA(LIST)(range(len(V)))
bases = [VV,EV,FV,CV]
orientedBoundary = signedCellularBoundaryCells(V,bases)
cells = [FV[f] if sign==1 else REVERSE(FV[f])  for (sign,f) in zip(*orientedBoundary)]
VIEW(EXPLODE(1.25,1.25,1.25)(MKPOLS((V,cells))))
