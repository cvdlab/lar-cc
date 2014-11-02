
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from bool1 import *

""" Definition of Boolean arguments """
n = 8
mod_1 = AA(LIST)(range(n)), [[2*k,2*k+1] for k in range(n/2)]
squares1 = INSR(larModelProduct)([mod_1,mod_1,mod_1])
mod_2 = AA(LIST)([0.5+k*2 for k in range(n/2)]),[[2*k,2*k+1] for k in range(n/4)]
squares2 = INSR(larModelProduct)([mod_2,mod_2,mod_2])

V1 = squares1[0]
V2 = squares2[0]
VV1 = AA(LIST)(range(len(V1)))
VV2 = AA(LIST)(range(len(V2)))
FV1 = larConvexFacets (*squares1)
FV2 = larConvexFacets (*squares2)
CV1 = squares1[1]
CV2 = squares2[1]

arg1 = V1,(VV1,FV1,CV1)
arg2 = V2,(VV2,FV2,CV2)

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

