""" A mesh model and various incidence operators """
import sys
sys.path.insert(0, 'lib/py/')
from larcc import *
from largrid import *

shape = [2,2,2]
V,(VV,EV,FV,CV) = larCuboids(shape,True)
"""
CV = [cell for cell in cellComplex if len(cell)==8]
FV = [cell for cell in cellComplex if len(cell)==4]
EV = [cell for cell in cellComplex if len(cell)==2]
VV = [cell for cell in cellComplex if len(cell)==1]
"""
VIEW(modelIndexing(shape))

CF = larCellFace(CV,FV)
CE = larCellFace(CV,EV)
FE = larCellFace(FV,EV)
