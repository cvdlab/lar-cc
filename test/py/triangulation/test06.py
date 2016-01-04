""" Generating the Triangulation of a set of non-intersecting cycles """
from larlib import *

sys.path.insert(0, 'test/py/triangulation/')
from test03 import *

triangleSet = larTriangulation( (V,EV) )

VIEW(STRUCT(AA(JOIN)(AA(AA(MK))(CAT(triangleSet)))))
VIEW(SKEL_1(STRUCT(AA(JOIN)(AA(AA(MK))(CAT(triangleSet))))))

"""
model = V,EV
W,FW = lar2boundaryPolygons(model)
polygons = [[W[u] for u in poly] for poly in FW]
VIEW(STRUCT(AA(POLYLINE)(polygons)))

triangleSet,triangledFace = [],[]
for polygon in polygons:  
    triangledPolygon = []
    polyline = []
    for p in polygon:
        polyline.append(Point(p[0],p[1]))
    cdt = CDT(polyline)

    triangles = cdt.triangulate()
    trias = [ [[t.a.x,t.a.y,0],[t.c.x,t.c.y,0],[t.b.x,t.b.y,0]] for t in triangles ]
    triangleSet += [AA(REVERSE)(trias)]
"""
