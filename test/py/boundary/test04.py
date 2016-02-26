""" 3D non-convex LAR cells """
from larlib import *
""" Input of a cellular 3-complex """
V,[VV,EV,FV,CV] = larCuboids([2,1,1],True)
struct = Struct([(V,FV,EV),t(.25,.25,0),s(.25,.5,2),(V,FV,EV)])

V,FV,EV = struct2Marshal(struct)
CF = AA(sorted)([[20,12,21,5,19,6],[27,1,5,28,13,23],[12,14,25,17,10,4],
[1,7,17,24,11,18],[30,29,26,16,8,22,10,11,4,18,24,25],[2,3,8,9,0,15]])
CV = [list(set(CAT([FV[f]  for f in faces]))) for faces in CF]

VV = AA(LIST)(range(len(V)))
hpc = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV,CV],hpc,0.6))

""" Visualization of the boundary 2-chain of a 3-complex """

V,BF,BE = larBoundary3(V,CV,FV,EV)(len(CV)*[1])
VIEW(STRUCT(MKTRIANGLES((V,BF,EV),color=True)))
VIEW(SKEL_1(STRUCT(MKTRIANGLES((V,BF,EV)) )))

boundaryEdges = chain2BoundaryChain(boundary2(FV,EV,VV))
edgeChain = boundaryEdges(29*[0]+[1]+[1]) 
VIEW(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,FV[29:31],[EV[e] for e in edgeChain]),color=True)))
VIEW(SKEL_1(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,FV[29:31],[EV[e] for e in edgeChain])))))

""" Visualization of a 3-chain of a 3-complex """

V,BF,BE = larBoundary3(V,CV,FV,EV)([0,0,0,0,1,1])
VIEW(STRUCT(MKTRIANGLES((V,BF,BE))))
VIEW(SKEL_1(STRUCT(MKTRIANGLES((V,BF,BE)) )))

