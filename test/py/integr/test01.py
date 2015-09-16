""" Integrals on the standard triangle """
from larlib import *

V = [[0,0,0],[1,0,0],[0,1,0]]
FV = [[0,1,2]]
P = (V,FV)
print II(P, 0, 0, 0, False)
