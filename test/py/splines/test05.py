""" Example of bilinear tensor product surface patch """
from larlib import *

controlpoints = [
   [[0,0,0],[2,-4,2]],
   [[0,3,1],[4,0,0]]]
dom = larDomain([20])
dom2D = larExtrude1(dom, 20*[1./20])
mapping = larBilinearSurface(controlpoints)
patch = larMap(mapping)(dom2D)
VIEW(STRUCT(MKPOLS(patch)))
