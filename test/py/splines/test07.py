""" Example of bilinear tensor product surface patch """
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/') 
from splines import *

controlpoints=[
   [[ 0,0,0],[0 ,3  ,4],[0,6,3],[0,10,0]],
   [[ 3,0,2],[2 ,2.5,5],[3,6,5],[4,8,2]],
   [[ 6,0,2],[8 ,3 , 5],[7,6,4.5],[6,10,2.5]],
   [[10,0,0],[11,3  ,4],[11,6,3],[10,9,0]]]
dom = larDomain([20])
dom2D = larExtrude1(dom, 20*[1./20])
mapping = larBicubicSurface(controlpoints)
patch = larMap(mapping)(dom2D)
VIEW(STRUCT(MKPOLS(patch)))
