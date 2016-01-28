
import sys
from larlib import *

sys.path.insert(0, 'test/py/bool/')
from test06 import *

""" From triples of points to LAR model """
FE = crossRelation(len(W),FW,EW)
triangleSet = boundaryTriangulation(W,FW,EW,FE)
TW,FT = triangleIndices(triangleSet,W)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,CAT(TW)))))
