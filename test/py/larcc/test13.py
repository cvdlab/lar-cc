""" Example of incidence chain computation """
import sys; sys.path.insert(0, 'lib/py/')
from larcc import *
from largrid import *

shape = (1,1,2) 
print "\n\nFor a better example provide a greater shape!"
V,bases = larCuboids(shape,True)

VV,EV,FV,CV = bases
incidence = incidenceChain([VV,EV,FV,CV])
relations = ["CF","FE","EV"]
for k in range(3):
   print "\n\n incidence", relations[k], "=\n", incidence[k],
print "\n\n"

submodel = SKEL_1(STRUCT(MKPOLS((V,EV))))
VIEW(larModelNumbering(V,[VV,EV,FV,CV],submodel,1))
