import sys; sys.path.insert(0, 'lib/py/')
from simplexn import *
from larcc import *

V,CV = larSimplexGrid1([10,10,3])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))
SK2 = (V,larSimplexFacets(CV))
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK2)))
_,FV = SK2
SK1 = (V,larSimplexFacets(FV))
_,EV = SK1
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK1)))

boundaryCells_2 = boundaryCells(CV,FV)
boundary = (V,[FV[k] for k in boundaryCells_2])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(boundary)))
print "\nboundaryCells_2 =\n", boundaryCells_2

boundaryCells_2 = signedBoundaryCells(V,CV,FV)
boundaryFV = [FV[-k] if k<0 else swap(FV[k]) for k in boundaryCells_2]

VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,boundaryFV))))
print "\nboundaryCells_2 =\n", boundaryFV

