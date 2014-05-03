""" Union of 2D non-structured grids """
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

import myfont
from myfont import *

import mapper
from mapper import *


import boolean2
from boolean2 import *

import lar2psm
from lar2psm import *

import myfont
from myfont import *

model1 = randomTriangulation(1000,2,'disk')
V1,CV1 = model1
VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(model1)+cellNames(model1,CV1,MAGENTA)))
model2 = randomTriangulation(1000,2,'cuboid')
V2,CV2 = model2
V2 = larScale( [2,2])(V2)
model2 = V2,CV2 
VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(model2)+cellNames(model2,CV2,RED)))
V, n1,n2,n12, B1,B2 = boolOps(model1,model2)
""" Visualization of first Boolean step  """
if n12==0:
   hpc0 = STRUCT([ COLOR(RED)(EXPLODE(1.5,1.5,1)(AA(MK)(V[:n1-n12]) )), 
            COLOR(CYAN)(EXPLODE(1.5,1.5,1)(AA(MK)(V[n1:]) )) ])
else:
   hpc0 = STRUCT([ COLOR(RED)(EXPLODE(1.5,1.5,1)(AA(MK)(V[:n1-n12]) )), 
            COLOR(CYAN)(EXPLODE(1.5,1.5,1)(AA(MK)(V[n1:]) )), 
            COLOR(WHITE)(EXPLODE(1.5,1.5,1)(AA(MK)(V[n1-n12:n1]) )) ])

# hpc1 = COLOR(RED)(EXPLODE(1.5,1.5,1)(MKPOLS((V,CV_un)) ))
# hpc2 = COLOR(CYAN)(EXPLODE(1.5,1.5,1)(MKPOLS((V,CV_int)) ))
# VIEW(STRUCT([hpc0, hpc1, hpc2]))

