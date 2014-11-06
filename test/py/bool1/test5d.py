
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from bool1 import *

V1 = [[0,0],[5,0],[5,5],[0,5]]
FV1 = [range(4)]
EV1 = [[0,1],[1,2],[2,3],[0,3]]
VV1 = AA(LIST)(range(len(V1)))

V2 = [[5,6],[10,6],[10,11],[5,11]]
FV2 = [range(4)]
EV2 = [[0,1],[1,2],[2,3],[0,3]]
VV2 = AA(LIST)(range(len(V2)))

arg1 = V1,(VV1,EV1,FV1)
arg2 = V2,(VV2,EV2,FV2)

""" Debug via visualization """
boolean = larBool(arg1,arg2)  

W,CW,chain,CX,FX,orientedBoundary = boolean("xor")
glass = MATERIAL([1,0,0,0.6,  0,1,0,0.6,  0,0,1,0.6, 0,0,0,0.6, 100])
VIEW(glass(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,chain)))))

if DEBUG:
   VIEW(SKEL_1(EXPLODE(1.1,1.1,1.1)(MKPOLS((W,orientedBoundary)))))
   
   W,CW,chain,CX,FX,orientedBoundary = boolean("union")
   VIEW(EXPLODE(1.1,1.1,1)(MKPOLS((W,chain))))
   VIEW(SKEL_1(EXPLODE(1.1,1.1,1.1)(MKPOLS((W,orientedBoundary)))))
   
   W,CW,chain,CX,FX,orientedBoundary = boolean("intersection")
   if chain != []:
      VIEW(EXPLODE(1.1,1.1,1)(MKPOLS((W,chain))))
      VIEW(SKEL_1(EXPLODE(1.1,1.1,1.1)(MKPOLS((W,orientedBoundary)))))
   
   W,CW,chain,CX,FX,orientedBoundary = boolean("difference")
   if chain != []:
      VIEW(EXPLODE(1.1,1.1,1)(MKPOLS((W,chain))))
      VIEW(SKEL_1(EXPLODE(1.1,1.1,1.1)(MKPOLS((W,orientedBoundary)))))
   
      VIEW(EXPLODE(1.1,1.1,1.1)(MKPOLS((W,CX))))

submodel = SKEL_1(STRUCT(MKPOLS((W,FX))))
VV = AA(LIST)(range(len(W)))
VIEW(larModelNumbering(1,1,1)(W,[VV,FX,CX],submodel,1))

