
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

