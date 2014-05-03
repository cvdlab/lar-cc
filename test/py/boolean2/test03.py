""" test example with general LAR cells for the boolean module """
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

V1 = [[0,5],[6,5],[0,2],[3,2],[6,2],[0,0],[3,0],[3,-2],[6,-2]]
CV1 = [[2,3,5,6],[0,1,2,3,4],[3,4,6,7,8]]
blue = V1,CV1
V2 = [[3,6],[7,6],[0,5],[3,5],[3,4],[7,4],[3,2],[7,2],[0,0],[3,0],[6,0],[6,2]]
CV2 = [[0,1,3,4,5],[2,3,4,6,8,9],[6,9,10,11],[4,5,6,7,11]]
red = V2,CV2
V, n1,n2,n12, B1,B2 = boolOps(red,blue)
CV = Delaunay(array(V)).vertices
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

