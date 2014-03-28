""" test program for the boolean module """
from pyplasm import *
from scipy import *
import os,sys

""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')

import boolean
from boolean import *

import lar2psm
from lar2psm import *

import simplexn
from simplexn import *

import larcc
from larcc import *

model1 = randomTriangulation(100,2,'disk')
VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(model1)))
model2 = randomTriangulation(100,2,'cuboid')
V2,CV2 = model2
V2 = scalePoints(V2, [2,2])
model2 = V2,CV2 
VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(model2)))
V,CV_un, CV_int, n1,n2,n12 = boolOps(model1,model2)
model = V,CV_int

hpc0 = STRUCT([ COLOR(RED)(EXPLODE(1.5,1.5,1)(AA(MK)(V[:n1-n12]) )), 
            COLOR(CYAN)(EXPLODE(1.5,1.5,1)(AA(MK)(V[n1:]) ))#, 
            #COLOR(WHITE)(EXPLODE(1.5,1.5,1)(AA(MK)(V[n1-n12:n1]) )) 
            ])
hpc1 = COLOR(RED)(EXPLODE(1.5,1.5,1)(MKPOLS((V,CV_un)) ))
hpc2 = COLOR(CYAN)(EXPLODE(1.5,1.5,1)(MKPOLS((V,CV_int)) ))
VIEW(STRUCT([hpc0, hpc1, hpc2]))
