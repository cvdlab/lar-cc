""" comparing oriented boundary and unoriented boundary extraction on a simple example """
import sys; sys.path.insert(0, 'lib/py/')
from largrid import *
from larcc import *

V,CV = larSimplexGrid1([1,1,1])
FV = larSimplexFacets(CV)

orientedBoundary = signedBoundaryCells(V,CV,FV)
orientedBoundaryFV = [FV[-k] if k<0 else swap(FV[k]) for k in orientedBoundary]
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,orientedBoundaryFV))))

BF = boundaryCells(CV,FV)
boundaryCellsFV = [FV[k] for k in BF]
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,boundaryCellsFV))))
