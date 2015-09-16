""" Boundary orientation of a random 2D cubical complex """
from larlib import *
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
orientedBoundary = signedCellularBoundaryCells(V,[VV,EV,FV])
cells = [EV[e] if sign==1 else REVERSE(EV[e]) for (sign,e) in zip(*orientedBoundary)]

# test model visualization
VIEW(STRUCT(MKPOLS((V,FV))))
VIEW(STRUCT(MKPOLS((V,EV))))
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,cells))))
VIEW(STRUCT(MKPOLS((V,cells))))
VIEW(mkSignedEdges((V,cells),2))
