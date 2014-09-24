
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from bool1 import *


n = 48
V1 = [[5*cos(angle*2*PI/n)+2.5, 5*sin(angle*2*PI/n)+2.5] for angle in range(n)]
FV1 = [range(n)]
EV1 = TRANS([range(n),range(1,n+1)]); EV1[-1] = [0,n-1]
VV1 = AA(LIST)(range(len(V1)))

V2 = [[4*cos(angle*2*PI/n), 4*sin(angle*2*PI/n)] for angle in range(n)]
FV2 = [range(n)]
EV2 = EV1
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

V,CV,vertDict,n1,n12,n2,BC,nbc1,nbc2 = makeCDC(arg1,arg2)      #<<<<<<<<<<<<<<<<

W,CW,VC,BCellCovering,cellCuts,boundary1,boundary2,BCW = makeSCDC(V,CV,BC,nbc1,nbc2)
assert len(VC) == len(V) 
assert len(BCellCovering) == len(BC)

submodel = STRUCT([ SKEL_1(STRUCT(MKPOLS((V,CV)))), COLOR(RED)(STRUCT(MKPOLS((V,BC)))) ])
dim = len(V[0])
if DEBUG: VIEW(STRUCT([ submodel,larModelNumbering(V,[VV,BC,CV],submodel,3)]))
if DEBUG: VIEW(EXPLODE(2,2,2)(MKPOLS((W,CW))))


if DEBUG: VIEW(STRUCT([ COLOR(GREEN)(EXPLODE(1.2,1.2,1)(MKPOLS((W,CAT(boundary1.values()))))), 
            COLOR(YELLOW)(EXPLODE(1.2,1.2,1)(MKPOLS((W,CAT(boundary2.values()))))) ]))


WW = AA(LIST)(range(len(W)))
FW = larConvexFacets (W,CW)
if len(CW)==4: FW=[[0,1],[1,2],[0,2],[0,3],[2,3],[2,4],[2,5],[3,4],[4,5]] #test5.py
_,EW = larFacets((W,FW), dim=2)

FWdict = dict()
for k,facet in enumerate (FW): FWdict[str(facet)] = k
for key,value in boundary1.items():
   value = [FWdict[str(facet)] for facet in value]
   boundary1[key] = value
for key,value in boundary2.items():
   value = [FWdict[str(facet)] for facet in value]
   boundary2[key] = value

glass = MATERIAL([1,0,0,0.2,  0,1,0,0.2,  0,0,1,0.1, 0,0,0,0.1, 100])
submodel = glass(STRUCT(MKPOLS((W,BCW))))
if DEBUG: VIEW(STRUCT([ submodel, larModelNumbering(W,[WW,FW,CW],submodel,1.5) ]))

if dim == 3: 
   _,EW = larFacets((W,FW), dim=2)
   bases = [WW,EW,FW,CW]
elif dim == 2: bases = [WW,FW,CW]
else: print "\nerror: not implemented\n"

coBoundaryMat = signedCellularBoundary(W,bases).T
boundaryMat = coBoundaryMat.T

CWbits = [[-1,-1] for k in range(len(CW))]
CWbits = cellTagging(boundary1,boundaryMat,CW,FW,W,BCW,CWbits,0)
CWbits = cellTagging(boundary2,boundaryMat,CW,FW,W,BCW,CWbits,1)

if DEBUG: VIEW(STRUCT(MKPOLS((W,BCW))))

"""
for k in range(1,len(CW)+1):
   VIEW(STRUCT([ STRUCT(MKPOLS((W,CW[:k]))), submodel,larModelNumbering(W,[BCW],submodel,3) ]))
"""

for cell in range(len(CW)):
   if CWbits[cell][0] == 1:
      CWbits = booleanChainTraverse(0,cell,W,CW,CWbits,1)      
   if CWbits[cell][0] == 0:
      CWbits = booleanChainTraverse(0,cell,W,CW,CWbits,0)
   if CWbits[cell][1] == 1:
      CWbits = booleanChainTraverse(1,cell,W,CW,CWbits,1)
   if CWbits[cell][1] == 0:
      CWbits = booleanChainTraverse(1,cell,W,CW,CWbits,0)

chain1,chain2 = TRANS(CWbits)
#VIEW(STRUCT([ submodel, STRUCT(MKPOLS((W,[cell for cell,c in zip(CW,chain1) if c==0] ))) ]) )
#VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((W,[cell for cell,c in zip(CW,chain2) if c==1] ))))

#VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((W,[cell for cell,c1,c2 in zip(CW,chain1,chain2) if c1*c2==1] ))))
VIEW(EXPLODE(1.1,1.1,1)(MKPOLS((W,[cell for cell,c1,c2 in zip(CW,chain1,chain2) if c1+c2==1] ))))
#VIEW(EXPLODE(1.1,1.1,1.1)(MKPOLS((W,[cell for cell,c1,c2 in zip(CW,chain1,chain2) if c1+c2>=1] ))))

submodel = STRUCT(MKPOLS((W,FW)))
VIEW(STRUCT([ submodel,larModelNumbering(W,[WW,FW,CW],submodel,2) ]))

W,CX = gatherPolytopes(W,CW,FW,boundaryMat,boundary1,boundary2)
print "\n CX =",CX
FX = larConvexFacets (W,CX)
VIEW(SKEL_1(EXPLODE(1.2,1.2,1)( MKPOLS((W,CX)) )))
VIEW(SKEL_1(EXPLODE(1.2,1.2,1)( MKPOLS((W,FX)) )))                      

