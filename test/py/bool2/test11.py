
""" Visualization of indices of the boundary triangulation """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *
sys.path.insert(0, 'test/py/bool2/')
from test06 import *

model = W,FW,EW
EF_angle = faceSlopeOrdering(model)

WW = AA(LIST)(range(len(W)))
submodel = SKEL_1(STRUCT(MKPOLS((W,EW))))
VIEW(larModelNumbering(1,1,1)(W,[WW,EW,FW],submodel,0.6))
