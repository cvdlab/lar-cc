""" Boundary of a 3-complex """
from larlib import *

V,[VV,EV,FV,CV] = larCuboids([1,1,1],True)
cube = Struct([ (V,FV,EV) ])
assembly = Struct([ cube, Struct([t(0,.5,0), r(PI/4,0,0), s(.5,.5,.5),cube]) ])

V,FV,EV = struct2Marshal(assembly)
VV = AA(LIST)(range(len(V)))
hpc = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],hpc,0.6))

CF = [[1,2,3,4,6,7],[0,1,2,3,4,5,6,7,8,9,10,11]]
CV = [list(set(CAT([FV[f]  for f in faces]))) for faces in CF]

V,BF,BE = larBoundary3(V,CV,FV,EV)([0,1])
VIEW(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,BF,BE)))) # ERROR in MKTRIANGLES with non-manifold face
VIEW(EXPLODE(1.2,1.2,1.2)(MKFACES((V,BF,EV))))
VIEW(SKEL_1(EXPLODE(1.2,1.2,1.2)(MKFACES((V,BF,EV)))))
