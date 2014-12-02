
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from bool1 import *

V1 = [[3,0],[11,0],[13,10],[10,11],[8,11],[6,11],[4,11],[1,10],[4,3],[6,4],
      [8,4],[10,3]]
FV1 = [[0,1,8,9,10,11],[1,2,11],[3,10,11],[4,5,9,10],[6,8,9],[0,7,8]]
EV1 = [[0,1],[0,7],[0,8],[1,2],[1,11],[2,11],[3,10],[3,11],[4,5],[4,10],[5,
      9],[6,8],[6,9],[7,8],[8,9],[9,10],[10,11]]
VV1 = AA(LIST)(range(len(V1)))

V2 = [[0,3],[14,2],[14,5],[14,7],[14,11],[0,8],[3,7],[3,5]]
FV2 = [[0,5,6,7],[0,1,7],[4,5,6],[2,3,6,7]]
EV2 = [[0,1],[0,5],[0,7],[1,7],[2,3],[2,7],[3,6],[4,5],[4,6],[5,6],[6,7]]
VV2 = AA(LIST)(range(len(V2)))

arg1 = V1,(VV1,EV1,FV1)
arg2 = V2,(VV2,EV2,FV2)

""" Debug via visualization """

V1,(VV1,EV1,FV1) = arg1
V2,(VV2,EV2,FV2) = arg2
glass = MATERIAL([1,0,0,0.3,  0,1,0,0.3,  0,0,1,0.3, 0,0,0,0.3, 100])

VIEW(STRUCT([
   glass(EXPLODE(1.1,1.1,1.1)(MKPOLS((V1,FV1)))),
   glass(EXPLODE(1.1,1.1,1.1)(MKPOLS((V2,FV2))))
]))

glass = MATERIAL([1,0,0,0.6,  0,1,0,0.6,  0,0,1,0.6, 0,0,0,0.6, 100])

boolean = larBool(arg1,arg2)  

W,CW,chain,CX,FX,orientedBoundary = boolean("xor")
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

