""" Boundary of a 2D cuboidal grid """
import sys;sys.path.insert(0, 'lib/py/')
from larcc import *

V,bases = larCuboids([6,6],True)
[VV,EV,FV] = bases
submodel = mkSignedEdges((V,EV))
VIEW(submodel)
VIEW(larModelNumbering(V,bases,submodel,1))

orientedBoundary = signedCellularBoundaryCells(V,bases)
submodel = mkSignedEdges((V,orientedBoundary))
VIEW(submodel)
VIEW(larModelNumbering(V,bases,submodel,1))
