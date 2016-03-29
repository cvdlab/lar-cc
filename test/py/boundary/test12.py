""" testing boundary operators (correct result) """
from larlib import *

V,[VV,EV,FV,CV] = larCuboids([1,1,1],True)
cell = (V,FV,EV)
cubeGrid =  Struct([Struct([ s(10,10,1),cell ])] + 3*[t(0,2,0), Struct(3*[ t(1.5,0,0), cell] )] ,"cubeGrid")
VIEW(STRUCT(MKPOLS(struct2lar(cubeGrid))))

V,FV,EV = struct2Marshal(cubeGrid)
csrmat,CF,faceCounter = larSignedBoundary3((V,FV,EV))
print csrmat.todense()
