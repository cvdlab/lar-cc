""" Boundary extraction of a portion of hollow sphere """
from pyplasm import *
from scipy import *
import os,sys

""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
import lar2psm
from lar2psm import *

import simplexn
from simplexn import *

import larcc
from larcc import *

import largrid
from largrid import *

import boolean2
from boolean2 import *


from mapper import *
from boolean import boolOps
model = larHollowSphere(0.8,1,PI/6,PI/4)([6,12,2])
V,FV = larHollowSphereFacets(0.8,1,PI/6,PI/4)([6,12,2])
print "\nV,FV =",(V,FV)
V,CV = model
print "\nV,CV =",(V,CV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))
BF = boundaryCells(CV,FV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,[FV[f] for f in BF]))))
