""" Generating the Triangulation of a set of non-intersecting cycles """
from larlib import *

sys.path.insert(0, 'test/py/inters/')
from test17 import *

model = V,EV
W,FW = lar2boundaryPolygons(model)
polygons = [[W[u] for u in poly] for poly in FW]
VIEW(STRUCT(AA(POLYLINE)(polygons)))

triangleSet = []  
for polygon in polygons:  
    polyline = [Point(p[0],p[1]) for p in polygon]
    cdt = CDT(polyline)
    triangles = cdt.triangulate()
    trias = [ [[t.c.x,t.c.y,0,1],[t.b.x,t.b.y,0,1],[t.a.x,t.a.y,0,1]] for t in triangles ]

    triangledFace += [[v[:-1] for v in triangle] for triangle in trias]
    triangleSet += [triangledFace]
