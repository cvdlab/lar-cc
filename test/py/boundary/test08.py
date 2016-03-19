""" Boundary of a 3-complex """
from larlib import *

V,[VV,EV,FV,CV] = larCuboids([1,1,1],True)
cube = Struct([ (V,FV,EV) ])
hole = Struct([t(0,.5,0), r(PI/4,0,0), s(1,.5/SQRT(2),.5/SQRT(2)),cube])
assembly = Struct([ cube, hole ])
assembly2 = Struct([ assembly, t(0,0,.5), s(0.5,1,1), hole ])

V,FV,EV = struct2Marshal(assembly2) # WRONG:  TODO: check ...
VV = AA(LIST)(range(len(V)))
hpc = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],hpc,0.7))

CF = [[1,3,6,14,17,18,19, 0,4,9,11,15,16, 2,5,7,8,10,13],[0,4,9,11,15,16],[2,5,7,8,10,13]]
CV = [list(set(CAT([FV[f]  for f in faces]))) for faces in CF]

V,BF,BE = larBoundary3(V,CV,FV,EV)([1,0,0])
VIEW(STRUCT(MKTRIANGLES((V,BF,BE),color=True))) 
VIEW(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,BF,BE),color=True))) 
VIEW(SKEL_1(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,BF,BE)))))
