
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
triangleSets = boundaryTriangulation(V,FV)
VIEW(EXPLODE(1.2,1.2,1.2)([STRUCT([MKPOL([tria,[[1,2,3,4]],None]) for tria in triangleSet]) for triangleSet in triangleSets]))

CF = AA(list)(CF)
CE = AA(list)(CE)

VIEW(EXPLODE(2,2,2) (AA(STRUCT)(AA(MKPOLS)( DISTL([V,[[EV[c] for c in cell] for cell in CE[:-1] ]])))  ))
VIEW(EXPLODE(2,2,2) (AA(STRUCT)(AA(MKPOLS)( DISTL([V,[[FV[c] for c in cell] for cell in CF ]])))  ))

models = DISTL([V,[[FV[c] for c in cell] for cell in CF ]])
models = [boundaryTriangulation(*model) for model in models]

def MKCELL(model): 
   return STRUCT([ STRUCT([MKPOL([tria,[[1,2,3,4]],None]) for tria in triangleSet]) 
          for triangleSet in model ])

VIEW(EXPLODE(1.5,1.5,1.5)(AA(MKCELL)([models[0],models[1],models[2]])))

WW = AA(LIST)(range(len(W)))
submodel = SKEL_1(STRUCT(MKPOLS((W,EW))))
VIEW(larModelNumbering(1,1,1)(W,[WW,EW,FW],submodel,0.6))
