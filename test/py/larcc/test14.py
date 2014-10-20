""" Boundary of a 2D cuboidal grid """
import sys;sys.path.insert(0, 'lib/py/')
from larcc import *

V,bases = larCuboids([6,6],True)
[VV,EV,FV] = bases
submodel = mkSignedEdges((V,EV))
VIEW(submodel)
VIEW(larModelNumbering(scalx=1,scaly=1,scalz=1)(V,bases,submodel,1))

orientedBoundary = signedCellularBoundaryCells(V,bases)
FV = [EV[e] if sign==1 else REVERSE(EV[e])  for (sign,e) in zip(*orientedBoundary)]
submodel = mkSignedEdges((V,FV))
VIEW(submodel)
VIEW(larModelNumbering(scalx=1,scaly=1,scalz=1)(V,bases,submodel,1))
