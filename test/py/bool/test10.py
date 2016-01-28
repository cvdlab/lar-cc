
""" Visualization of indices of the boundary triangulation """
from larlib import *

sys.path.insert(0, 'test/py/bool/')
from test09 import *

model = W,FW,EW
FE = crossRelation(len(W),FW,EW)
EF_angle, ET,TV = faceSlopeOrdering(model,FE)

WW = AA(LIST)(range(len(W)))
submodel = SKEL_1(STRUCT(MKPOLS((W,CAT(TW)))))
VIEW(larModelNumbering(1,1,1)(W,[WW,EW,CAT(TW)],submodel,0.6))
