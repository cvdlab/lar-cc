""" Correction to boundary operator """
from larlib import *

V,[VV,EV,FV] = larCuboids([2,1],True)
complex = Struct([(V,FV,EV), t(.5,.25), s(.5,.5), (V,FV,EV)])
V,FV,EV = struct2lar(complex)
lines = [[V[v] for v in edge] for edge in EV]
V,FV,EV,polygons = larFromLines(lines)

VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,.5)) 

csrBoundaryMat = boundary(FV,EV)
print "wrong boundary matrix =",csrBoundaryMat.todense()
csrBoundaryMat = larUnsignedBoundary2(FV,EV,VV)  # <<<<< NOTE !!
print "right boundary matrix =",csrBoundaryMat.todense()

""" View boundary chain """
def viewBoundaryChain(larModel):
    V,FV,EV = larModel
    VV = AA(LIST)(range(len(V)))
    def viewBoundaryChain0 (chain):
        BE = chain2BoundaryChain(larUnsignedBoundary2(FV,EV,VV))(chain)
        if chain == [1,1,1,1]: assert BE == [0,4,5,7,9,10]
        if chain == [0,0,1,1]: assert BE == [0,2,3,4,5,6,7,8,9,10,12,13]
        BEV = [EV[e] for e in BE]
        if chain == [1,1,1,1]: assert BEV == [(0,1),(5,6),(7,8),(5,7),(0,6),(1, 8)]
        submodel = STRUCT(MKPOLS((V,BEV)))
        return submodel
    return viewBoundaryChain0
V = [[0.0,1.0],[1.0,1.0],[1.0,0.75],[1.0,0.25],[0.5,0.25],[1.0,0.0],
[0.0,0.0],[2.0,0.0],[2.0,1.0],[0.5,0.75],[1.5,0.25],[1.5,0.75]]

FV = [[3,2,11,10],[9,2,3,4],[1,2,9,4,3,5,6,0],[1,8,7,5,3,10,11,2]]

EV = [(0,1),(1,2),(4,9),(2,9),(5,6),(7,8),(3,10),(5,7),(2,11),(0,6),
(1,8),(2,3),(10,11),(3,4),(3,5)]


submodel = viewBoundaryChain((V,FV,EV))([1,1,1,1])
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,.5)) 
submodel = viewBoundaryChain((V,FV,EV))([0,0,1,1])
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,.5)) 
