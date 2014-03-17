import sys,os
import scipy.misc, numpy, pickle
from numpy.random import randint
from pyplasm import *

""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')

import largrid
from largrid import *


import morph
from morph import *
 
shape = 64,64
structure = 8,8
assert len(shape) == len(structure)
imageVerts = larImageVerts(shape)
_, skeletons, operators = imageChainComplex (shape)
image_array = randomImage(shape, structure, 0.05)
minPoint, maxPoint = (0,0), (64,64)
window = minPoint, maxPoint
segmentChain = setMaskWindow(window,image_array)
   
solid = visImageChain (shape,segmentChain, imageVerts, skeletons)
b_rep,boundaryChain = imageChainBoundary(shape, operators)(2)(segmentChain)

stepwiseTest = True
if stepwiseTest:
   model = testAlgebraicMorphologyStepByStep(solid, b_rep, 
            boundaryChain, imageVerts, skeletons)
else:
   model = testAlgebraicMorphology(solid, b_rep, 
            boundaryChain, imageVerts, skeletons)
VIEW(STRUCT(MKPOLS(model)))

V = model[0]
M = AA(tuple)(model[1])
S = AA(tuple)(solid[1])
B = AA(tuple)(b_rep[1])
D = list(set(S).union(M))
E = list(set(S).difference(M))
M,S,B,D,E = (AA(AA(list)))([M,S,B,D,E])
M2 = STRUCT(MKPOLS((V,M)))
S2 = STRUCT(MKPOLS((V,S)))

S2 = COLOR(CYAN)(STRUCT(MKPOLS((V,S))))
B1 = COLOR(MAGENTA)(STRUCT(MKPOLS((V,B))))
D2 = COLOR(YELLOW)(STRUCT(MKPOLS((V,D))))
E2 = COLOR(WHITE)(STRUCT(MKPOLS((V,E))))
VIEW(STRUCT([D2,S2,E2,B1]))

VIEW(STRUCT([D2,S2,B1]))
VIEW(STRUCT([S2,E2,B1]))

