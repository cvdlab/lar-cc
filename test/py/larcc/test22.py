""" Testing correction to boundary operator for general (non-convex) LAR """
from larlib import *

""" Input and visualization of a general cellular complex """

mod1 = larQuote1D([0.2,-0.2,0.2])
mod2 = larQuote1D([0.5,0.5])
mod3 = larQuote1D([0.5])
mod = larModelProduct([larModelProduct([mod2,mod1]),mod3])

mx = larQuote1D([1,1])
my = larQuote1D([1])
m = larModelProduct([larModelProduct([mx,my]),my])

complex = Struct([m, t(.5,.2,.25), mod])
V,CV0 = struct2lar(complex)
CV = copy(CV0)
CV[1] += CV[0] + CV[4]
CV[3] += CV[2] + CV[5]
VIEW(SKEL_1(STRUCT( MKPOLS((V,CV)))))

""" Visualization of the skeletons """

def cubeFV(verts): 
    [a,b,c,d,  e,f,g,h] = verts
    return [[a,b,c,d],[e,f,g,h],[a,c,e,g],[b,d,f,h],[a,b,e,f],[c,d,g,h]]

FV = COMP([sorted,set,AA(tuple),CAT,AA(cubeFV)])(CV0)
V,EV = larFacets((V,FV),dim=2)
VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV,CV],submodel,.25)) 
FV[5] += FV[16] + FV[26]
assert 30==len(boundaryCells(CV0,FV))  # assembly of 3 unrelated 3-complexes

""" Computation of corrected boundary operator """

totalChain = len(CV)*[1]
BF = chain2BoundaryChain(larUnsignedBoundary2(CV,FV,EV))(totalChain)
BFV = [FV[f] for f in BF]
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,BFV))))
glass = MATERIAL([1,0,0,0.1,  0,1,0,0.1,  0,0,1,0.1, 0,0,0,0.1, 100])
VIEW(STRUCT(AA(glass)(MKPOLS((V,BFV)))))
VIEW(glass(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,BFV)))))
V,BEV = larFacets((V,BFV),dim=1)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,[e for e in BEV if len(e)==2]))))

""" Visualization of some boundary chains """

BF = chain2BoundaryChain(larUnsignedBoundary2(CV,FV,EV))([1,0,1,1,0,0])
BFV = [FV[f] for f in BF]
V,BEV = larFacets((V,BFV),dim=1)
BEV = [e for e in BEV if len(e)==2]
VIEW(glass(STRUCT(MKPOLS((V,BFV)))))
submodel = STRUCT(MKPOLS((V,BEV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,BEV,BFV],submodel,.25)) 


obj = Struct([(V,BFV,BEV)])
W,FW,EW = struct2lar(obj)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,EW))))
VIEW(SKEL_1(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((W,FW,EW)))))


