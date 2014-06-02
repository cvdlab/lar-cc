import sys; sys.path.insert(0, 'lib/py/')
from simplexn import *
from larcc import *
V,CV = larSimplexGrid1([4,4,4])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))

FV = larSimplexFacets(CV)
EV = larSimplexFacets(FV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))

csrSignedBoundaryMat = signedSimplicialBoundary (CV,FV)
boundaryCells_2 = signedBoundaryCells(V,CV,FV)
boundaryFV = [FV[-k] if k<0 else swap(FV[k]) for k in boundaryCells_2]
boundary = (V,boundaryFV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(boundary)))

