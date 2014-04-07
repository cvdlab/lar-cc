""" Circumference of unit radius """
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
"""
model = checkModel(larCircle(1)())
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(model)))
model = checkModel(larDisk(1)([36,4]))
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(model)))
model = checkModel(larRing([.9, 1.])([36,2]))
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(model)))
model = checkModel(larCylinder([.5,2.])([32,1]))
VIEW(STRUCT(MKPOLS(model)))
model = checkModel(larSphere(1)())
VIEW(STRUCT(MKPOLS(model)))
model = larBall(1)()
VIEW(STRUCT(MKPOLS(model)))
model = larRod([.25,2.])([32,1])
VIEW(STRUCT(MKPOLS(model)))
model = checkModel(larToroidal([0.5,1])())
VIEW(STRUCT(MKPOLS(model)))
model = checkModel(larCrown([0.125,1])([8,48]))
VIEW(STRUCT(MKPOLS(model)))
model = larPizza([0.05,1])([8,48])
VIEW(STRUCT(MKPOLS(model)))
"""
model = checkModel(larTorus([0.5,1])())
VIEW(STRUCT(MKPOLS(model)))
