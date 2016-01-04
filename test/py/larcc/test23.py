""" 2-cell with high topological genus """
from larlib import *

""" 2-cell with high topological genus 0. """
side = QUOTE(5*[1,-1])
holes = PROD([side,side])
V,FV,_ = UKPOL(holes)
FV = [[v-1 for v in f] for f in FV]
EV = CAT([[[v,f[(k+1)%4]] for k,v in enumerate(f+[f[0]][:-1])] for f in FV])
VIEW(STRUCT(MKPOLS((V,EV))))

(W,_) = larBox((-1,-1),(10,10))
(_,(_,EW,FW)) = larCuboids((1,1),True)
complex = Struct([(W,FW,EW),(V,FV,EV)])
V,FV,EV = struct2lar(complex)

""" 2-cell with high topological genus 1. """
FV = [sorted(CAT(FV))] 
# EV = sorted(AA(sorted)(EV))  # p2t bug if uncommented: CHECK
VIEW(STRUCT(MKPOLS((V,EV))))

VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,1.5)) 

""" 2-cell with high topological genus 2. """

# LAR data structures
csrBoundaryMat = boundary(FV,EV)
boundaryChain = chain2BoundaryChain(csrBoundaryMat)([1])
triangleSet = larTriangulation( (V,EV) )

# PyPLaSM data structures
hpcChain = AA(JOIN)(AA(AA(MK))(CAT(triangleSet)))
hpcChainBoundary = AA(COLOR(RED))(MKPOLS((V,[EV[e] for e in boundaryChain])))

VIEW(STRUCT( hpcChain + hpcChainBoundary ))
VIEW(EXPLODE(1.2,1.2,1.2)( hpcChain + hpcChainBoundary ))

