""" test program for the boolean module """
from pyplasm import *
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

import largrid
from largrid import *

blue = larSimplexGrid([2,4])
red = larSimplexGrid([4,3])
VIEW(STRUCT([
COLOR(RED)(EXPLODE(1.2,1.2,1)(MKPOLS(red))),
COLOR(BLUE)(EXPLODE(1.2,1.2,1)(MKPOLS(blue)))
]))

V,CV_un, CV_int, n1,n2,n12 = boolOps(red,blue)
CV = Delaunay(array(V)).vertices

if n12==0:
   hpc0 = STRUCT([ COLOR(RED)(EXPLODE(1.5,1.5,1)(AA(MK)(V[:n1-n12]) )), 
            COLOR(CYAN)(EXPLODE(1.5,1.5,1)(AA(MK)(V[n1:]) )) ])
else:
   hpc0 = STRUCT([ COLOR(RED)(EXPLODE(1.5,1.5,1)(AA(MK)(V[:n1-n12]) )), 
            COLOR(CYAN)(EXPLODE(1.5,1.5,1)(AA(MK)(V[n1:]) )), 
            COLOR(WHITE)(EXPLODE(1.5,1.5,1)(AA(MK)(V[n1-n12:n1]) )) ])

hpc1 = COLOR(RED)(EXPLODE(1.5,1.5,1)(MKPOLS((V,CV_un)) ))
hpc2 = COLOR(CYAN)(EXPLODE(1.5,1.5,1)(MKPOLS((V,CV_int)) ))
VIEW(STRUCT([hpc0, hpc1, hpc2]))
