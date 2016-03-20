""" Boundary of a 3-complex """
from larlib import *

V,[VV,EV,FV,CV] = larCuboids([1,1,1],True)
cube = Struct([ (V,FV,EV) ])
hole = Struct([t(0,.5,0), r(PI/4,0,0), s(.5,.5,.5),cube])
assembly = Struct([ cube, hole, t(0,0,SQRT(0.5)), hole ])

V,FV,EV = struct2Marshal(assembly) # WRONG:  TODO: check ...
VV = AA(LIST)(range(len(V)))
hpc = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[[],[],FV],hpc,0.6))

CF = [[4,5,7,16,17,19,20],[3,8,6,12,11,13],[0,1,10,20],[]]
CV = [list(set(CAT([FV[f]  for f in faces]))) for faces in CF]

V,BF,BE = larUnsignedBoundary3(V,CV,FV,EV)([0,1,1,0])
VIEW(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,BF,BE)))) # ERROR in MKTRIANGLES with non-manifold face
VIEW(EXPLODE(1.2,1.2,1.2)(MKFACES((V,BF,EV))))
VIEW(SKEL_1(EXPLODE(1.2,1.2,1.2)(MKFACES((V,BF,EV)))))
