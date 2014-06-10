""" Boundary orientation of a random 2D cubical complex """
import sys;sys.path.insert(0, 'lib/py/')
from scipy import linalg
from larcc import *
from random import random

# test model generation
shape = 20,20
V,FV = larCuboids(shape)
cellSpan = prod(shape)
fraction = 0.5
remove = [int(random()*cellSpan) for k in range(int(cellSpan*fraction)) ]
FV = [FV[k] for k in range(cellSpan) if not (k in remove)]
_,EV = larCuboidsFacets((V,FV))
VV = AA(LIST)(range(len(V)))
orientedBoundaryCells = signedCellularBoundaryCells(V,[VV,EV,FV])

# test model visualization
VIEW(STRUCT(MKPOLS((V,FV))))
VIEW(STRUCT(MKPOLS((V,EV))))
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,orientedBoundaryCells))))
VIEW(STRUCT(MKPOLS((V,orientedBoundaryCells))))
VIEW(mkSignedEdges((V,orientedBoundaryCells),2))
