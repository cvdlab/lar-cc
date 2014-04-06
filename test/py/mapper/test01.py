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


import mapper
from mapper import larCircle
model = larCircle(1)()
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
