""" A mesh model and various incidence operators """

from larlib import *

shape = [2,2,2]
V,(VV,EV,FV,CV) = larCuboids(shape,True)
VIEW(modelIndexing(shape))

CF = larCellFace(CV,FV)
CE = larCellFace(CV,EV)
FE = larCellFace(FV,EV)
