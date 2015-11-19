""" Generating the LAR of a set of non-intersecting cycles """
from larlib import *

sys.path.insert(0, 'test/py/inters/')
from test17 import *

W,EW = boundaryCycles2vertexPermutation( (V,EV) )
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,EW))))
