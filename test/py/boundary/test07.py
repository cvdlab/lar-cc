""" Boundary of a 3-complex """
from larlib import *

V,[VV,EV,FV,CV] = larCuboids([1,1,1],True)
cube = Struct([ (V,FV,EV) ])
hole = Struct([t(0,.5,0), r(PI/4,0,0), s(1,.5/SQRT(2),.5/SQRT(2)),cube])
assembly = Struct([ cube, hole ])

V,FV,EV = struct2Marshal(assembly) # WRONG:  TODO: check ...
VV = AA(LIST)(range(len(V)))
hpc = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],hpc,0.6))

CF = [[1,3,6,7,12,11],[0,2,4,5,9,8]]
CV = [list(set(CAT([FV[f]  for f in faces]))) for faces in CF]

V,BF,BE = larBoundary3(V,CV,FV,EV)([1,0])
VIEW(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,BF,EV),color=True)))
