""" Arrangements with non-contractible cells """
from larlib import *

V,[VV,EV,FV,CV] = larCuboids([1,1,1],True)
cube1 = Struct([(V,FV,EV)],"cube1")
cube2 = Struct([t(.25,.25,-1),s(.5,.5,3),(V,FV,EV)],"cube2")
V,FV,EV = struct2lar(Struct([cube1,cube2]))

VIEW(STRUCT(MKPOLS((V,FV,EV))))
V,CV,FV,EV,CF,CE,COE,FE = thePartition(V,FV,EV)
