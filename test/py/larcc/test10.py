""" A mesh model and various incidence operators """
import sys; sys.path.insert(0, 'lib/py/')
from larcc import *
from largrid import *

shape = [2,2,2]
V,(VV,EV,FV,CV) = larCuboids(shape,True)
VIEW(modelIndexing(shape))

CF = larCellFace(CV,FV)
CE = larCellFace(CV,EV)
FE = larCellFace(FV,EV)
