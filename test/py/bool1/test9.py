
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from bool1 import *

V1 = [[0,0],[15,0],[15,14],[0,14]]
FV1 = [range(4)]
EV1 = [[0,1],[1,2],[2,3],[0,3]]
VV1 = AA(LIST)(range(len(V1)))

V2 = [[1,1],[7,1],[7,6],[1,6], [8,1],[14,1],[14,7],[8,7], [1,7],[7,7],[7,13],[1,13], [8,8],[14,8],[14,13],[8,13]]
FV2 = [range(4),range(4,8),range(8,12),range(12,16)]
EV2 = [[0,1],[1,2],[2,3],[0,3], [4,5],[5,6],[6,7],[4,7], [8,9],[9,10],[10,11],[8,11], [12,13],[13,14],[14,15],[12,15]]
VV2 = AA(LIST)(range(len(V2)))

arg1 = V1,(VV1,EV1,FV1)
arg2 = V2,(VV2,EV2,FV2)

V1,basis1 = arg1
V2,basis2 = arg2
cells1 = basis1[-1]
cells2 = basis2[-1]

if DEBUG: VIEW(STRUCT(MKPOLS((V1,basis1[1])) + MKPOLS((V2,basis2[1]))))

model1,model2 = (V1,cells1),(V2,cells2)
V, CV1,CV2, n1,n12,n2 = mergeVertices(model1,model2)  #<<<<<<<<<<<<<<<<

submodel = SKEL_1(STRUCT(MKPOLS((V,CV1+CV2)))) 
VV = AA(LIST)(range(len(V)))
if DEBUG: VIEW(STRUCT([ submodel,larModelNumbering(V,[VV,_,CV1+CV2],submodel,3)]))

V,CV,vertDict,n1,n12,n2,BC = makeCDC(arg1,arg2)    #<<<<<<<<<<<<<<<<

W,CW,VC,BCellCovering,cellCuts,BCfrags = makeSCDC(V,CV,BC)

assert len(VC) == len(V) 
assert len(BCellCovering) == len(BC)

submodel = STRUCT([ SKEL_1(STRUCT(MKPOLS((V,CV)))), COLOR(RED)(STRUCT(MKPOLS((V,BC)))) ])
dim = len(V[0])
VIEW(STRUCT([ submodel,larModelNumbering(V,[VV,BC,CV],submodel,3)]))
VIEW(EXPLODE(2,2,2)(MKPOLS((W,CW))))
"""
for k in range(1,len(CW)+1):
   VIEW(STRUCT([ STRUCT(MKPOLS((W,CW[:k]))), submodel,larModelNumbering(V,[VV,BC,CV],submodel,3) ]))
"""


WW = AA(LIST)(range(len(W)))
FW = larConvexFacets (W,CW)
#submodel = SKEL_1(STRUCT(MKPOLS((W,CW))))
#VIEW(larModelNumbering(W,[WW,FW,CW],submodel,3))

VIEW(EXPLODE(1.5,1.5,1)(MKPOLS((W,FW))))



