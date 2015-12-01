""" Explicit LAR of 2-cell with high topological genus """
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


FV = [CAT(FV)] + sorted(AA(list)(FV))[1:]
# EV = sorted(AA(sorted)(EV))  # p2t bug if uncommented: CHECK
model = V,FV,EV

VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,1.5)) 

# LAR data structures
chain = [1]+[0]*(len(FV)-1)
outModel,triangleSet = larComplexChain(model)(chain)

viewLarComplexChain(model)([0]+[1]*(len(FV)-1))
viewLarComplexChain(model)([1]+[0]*(len(FV)-1))
viewLarComplexChain(model)([1]+[0]*5 + [1]*(len(FV)-6))
viewLarComplexChain(model)([1]*(len(FV)))
