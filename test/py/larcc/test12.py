""" comparing edge orientation and oriented boundary extraction """
import sys; sys.path.insert(0, 'lib/py/')
from largrid import *
from larcc import *

V,FV = larSimplexGrid1([5,5])
EV = larSimplexFacets(FV)
VIEW(mkSignedEdges((V,EV)))

orientedBoundary = signedBoundaryCells(V,FV,EV)
orientedBoundaryEV = [EV[-k] if k<0 else swap(EV[k]) for k in orientedBoundary]
VIEW(mkSignedEdges((V,orientedBoundaryEV)))
