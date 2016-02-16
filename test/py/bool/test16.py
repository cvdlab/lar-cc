""" Generation of signed boundary operator of a general LAR complex """
from larlib import *

V = [[0.25,0.25,-5e-05],[0.25,0.75,-5e-05],[0.75,0.75,-5e-05],[0.75,0.25,-5e-05],[
1.0,0.0,0.0],[0.0,0.0,0.0],[1.0,1.0,0.0],[0.0,1.0,0.0],[0.25,0.25,1.0],[0.25,0.25,
2.0],[0.25,0.75,2.0],[0.25,0.75,1.0],[0.25,0.75,-1.0],[0.25,0.25,-1.0],[0.75,0.75,
-1.0],[0.75,0.25,-1.0],[0.75,0.25,1.0],[0.75,0.75,1.0],[1.0,0.0,1.0],[0.0,0.0,
1.0],[1.0,1.0,1.0],[0.0,1.0,1.0],[0.75,0.75,2.0],[0.75,0.25,2.0]]

FV = [(2,3,16,17),(6,7,20,21),(12,13,14,15),(0,1,8,11),(1,2,11,17),(0,1,12,13),
(4,6,18,20),(5,7,19,21),(0,3,13,15),(0,3,8,16),(0,1,2,3),(10,11,17,22),(2,3,14,
15),(8,9,16,23),(8,11,16,17),(1,2,12,14),(16,17,22,23),(4,5,18,19),(8,9,10,11),(
9,10,22,23),(0,1,2,3,4,5,6,7),(8, 11,16,17,18,19,20,21)]

EV =[(3,15),(7,21),(10,11),(4,18),(12,13),(5,19),(8,9),(18,19),(22,23),(0,3),(1,11),
(16,17),(0,8),(6,7),(20,21),(3,16),(10,22),(18,20),(19,21),(1,2),(12,14),(4,5),(
8,11),(13,15),(16,23),(14,15),(11,17),(17,22),(2,14),(2,17),(0,1),(9,10),(8,16),
(4,6),(1,12),(5,7),(0,13),( 9,23),(6,20),(2,3)]

VV = AA(LIST)(range(len(V)))
submodel = SKEL_1(STRUCT(MKPOLS((V,EV))))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,0.6))

FE = crossRelation(FV,EV,VV)
EF_angle, ET,TV,FT = faceSlopeOrdering((V,FV,EV),FE)
EW = extendEV(EV,ET,TV)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,TV))))
            
submodel = SKEL_1(STRUCT(MKPOLS((V,EW))))
VIEW(larModelNumbering(1,1,1)(V,[VV,EW,TV],submodel,0.6))

triaModel, larModel = (TV,EW), (FV,EV)
op = larSignedBoundary(larModel,triaModel,FT)
print op.todense()
