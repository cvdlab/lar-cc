""" Example of nested structures with translation and rotations """
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
square = larCuboids([1,1])
square = Model(square,2)
table = larApply( t(-.5,-.5) )(square)
chair = Struct([ t(.75, 0), s(.35,.35), table ])
struct = Struct( [t(2,1)] + [table] + 4*[r(PI/2), chair])
scene = evalStruct(struct)
VIEW(SKEL_1(STRUCT(CAT(AA(MKPOLS)(scene)))))
