
""" Visualization of indices of the boundary triangulation """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *
sys.path.insert(0, 'test/py/bool2/')
from test06 import *

global count
count = 0

model = W,FW,EW
EF_angle = faceSlopeOrdering(model)

V,CV,FV,EV,CF,CE = facesFromComponents(model)
CF = AA(list)(CF)
CE = AA(list)(CE)
VIEW(EXPLODE(2,2,2)(MKPOLS((V,[CAT([FV[c] for c in cell]) for cell in CF]))))

VIEW(EXPLODE(2,2,2) (AA(STRUCT)(AA(MKPOLS)( DISTL([V,[[EV[c] for c in cell] for cell in [CE[-1]] ]])))  ))
VIEW(EXPLODE(2,2,2) (AA(STRUCT)(AA(MKPOLS)( DISTL([V,[[FV[c] for c in cell] for cell in CF]])))  ))

WW = AA(LIST)(range(len(W)))
submodel = SKEL_1(STRUCT(MKPOLS((W,EW))))
VIEW(larModelNumbering(1,1,1)(W,[WW,EW,FW],submodel,0.6))
