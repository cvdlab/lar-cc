""" 3D non-convex LAR cells """
from larlib import *

V,[VV,EV,FV,CV] = larCuboids([2,1,1],True)
mod1 = Struct([(V,FV,EV),t(.25,.25,0),s(.25,.5,2),(V,FV,EV)])
W,FW,EW = struct2lar(mod1)

quadArray = [[W[v] for v in face] for face in FW]
parts = boxBuckets3d(containmentBoxes(quadArray))
Z,FZ,EZ = spacePartition(W,FW,EW, parts)
Z,FZ,EZ = larSimplify((Z,FZ,EZ),radius=0.0001)
V,FV,EV = Z,FZ,EZ

CF = AA(sorted)([[20,12,21,5,19,6],[27,1,5,28,13,23],[12,14,25,17,10,4],
[1,7,17,24,11,18],[30,29,26,16,8,22,10,11,4,18,24,25],[2,3,8,9,0,15]])

CV = [list(set(CAT([FV[f]  for f in faces]))) for faces in CF]

VV = AA(LIST)(range(len(V)))
hpc = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV,CV],hpc,0.6))

BF = boundaryCells2(CV,FV,EV)
VIEW(SKEL_1(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,[FV[f] for f in BF],EV)))))
VIEW(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,[FV[f] for f in BF],EV))))

boundaryChain = chain2BoundaryChain(boundary1(FV,EV,VV))
faceChain = boundaryChain(29*[0]+[1]) + boundaryChain(30*[0]+[1])
VIEW(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES(
    (V,[FV[29],FV[30]],[EV[e] for e in faceChain]))))
VIEW(SKEL_1(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES(
    (V,[FV[29],FV[30]],[EV[e] for e in faceChain])))))

cellBoundary = chain2BoundaryChain(boundary2(CV,FV,EV))([0,0,0,0,1,0])
faceChain = set()
for f in cellBoundary: 
    faceCoords = len(FV)*[0]
    faceCoords[f] = 1
    faceChain = faceChain.union(boundaryChain(faceCoords))
VIEW(STRUCT(MKTRIANGLES(
    (V,[FV[f] for f in cellBoundary],[EV[e] for e in faceChain]))))
VIEW(SKEL_1(STRUCT(MKTRIANGLES(
    (V,[FV[f] for f in cellBoundary],[EV[e] for e in faceChain])))))
