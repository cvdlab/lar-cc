""" Integrals on the standard triangle """
import sys; sys.path.insert(0, 'lib/py/')
from integr import *

V = [[0,0,0],[1,0,0],[0,1,0]]
FV = [[0,1,2]]
P = (V,FV)
print II(P, 0, 0, 0)
