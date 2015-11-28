""" Correction to boundary operator """
from larlib import *

V,[VV,EV,FV] = larCuboids([2,1],True)
complex = Struct([(V,FV,EV), t(.5,.25), s(.5,.5), (V,FV,EV)])
V,FV,EV = struct2lar(complex)
lines = [[V[v] for v in edge] for edge in EV]
V,FV,EV = larFromLines(lines)

VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,.5)) 

csrBoundaryMat = boundary(FV,EV)
print "wrong boundary matrix =",csrBoundaryMat.todense()
csrBoundaryMat = boundary1(FV,EV,VV)  % <<<<< NOTE !!
print "right boundary matrix =",csrBoundaryMat.todense()
