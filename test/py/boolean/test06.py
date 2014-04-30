""" test program for the boolean module """
from pyplasm import *
from scipy import *
import os,sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from lar2psm import *
from simplexn import *
from larcc import *
from largrid import *
from myfont import *
from mapper import *

from mapper import *
from boolean import boolOps
blue = larHollowCyl(0.8,1,1,angle=PI/4)([6,2,5])
VIEW(STRUCT(MKPOLS(blue)))
V1,FV1 = larHollowCylFacets(0.8,1,1,angle=PI/4)([6,2,5])
assert blue[0]==V1
print "*** len(V1) =",len(V1)
red = larHollowSphere(0.8,1,PI/6,PI/4)([6,12,2])
VIEW(STRUCT(MKPOLS(red)))
V2,FV2= larHollowSphereFacets(0.8,1,PI/6,PI/4)([6,12,2])
assert red[0]==V2
print "*** len(V2) =",len(V2)
V, n1,n2,n12,BV1,BV2 = boolOps(blue,red,'cuboid',FV1,FV2)
