""" Integrals on the standard 3D cube """
import sys; sys.path.insert(0, 'lib/py/')
from integr import *

V = [[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]]

FV = [[1,0,2],[0,1,4],[2,0,4],[1,2,3],[1,3,5],[4,1,5],[3,2,6],[2,4,6],
     [5,3,7],[3,6,7],[4,5,6],[6,5,7]]

P = (V,FV)
print Volume(P)
print Centroid(P)

""" changing the boundary orientation changes the sign of volume,
    but not changes the centroid """

P = (V,AA(REVERSE)(FV))
print Volume(P)
print Centroid(P)
