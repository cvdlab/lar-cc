""" Boundary extraction of a portion of hollow sphere """

from larlib import *

# broken, TODO: debug

model = larHollowSphere([0.8,1])([6,12,1])
# V,FV = larHollowSphereFacets(0.8,1,PI/6,PI/4)([6,12,2])
print "\nV,FV =",(V,FV)
V,CV = model
print "\nV,CV =",(V,CV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))
BF = boundaryCells(CV,FV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,[FV[f] for f in BF]))))
