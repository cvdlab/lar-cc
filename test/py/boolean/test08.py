
import sys
from larlib import *

sys.path.insert(0, 'test/py/boolean/')
from test06 import *

""" From triples of points to LAR model """
WW = AA(LIST)(range(len(W)))
FE = crossRelation(FW,EW,WW)
triangleSet = boundaryTriangulation(W,FW,EW,FE)
TW,FT = triangleIndices(triangleSet,W)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,CAT(TW)))))
