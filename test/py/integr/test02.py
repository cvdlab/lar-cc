""" Integrals on the standard triangle """
from larlib import *

V = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
FV = [[1,2,3],[0,3,2],[0,1,3],[0,1,2]]
P = (V,FV)
print Volume(P)
