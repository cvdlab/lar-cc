""" 3D boundary triangulation of the space partition """
import sys
sys.path.insert(0, 'test/py/bool2/')
from test06 import *
from support import PolygonTessellator,vertex

def orientTriangle(pointTriple):
    v1 = array(pointTriple[1])-pointTriple[0]
    v2 = array(pointTriple[2])-pointTriple[0]
    if cross(v1,v2)[2] < 0: return REVERSE(pointTriple)
    else: return pointTriple

def boundaryTriangulation(W,FW):
    triangles = []
    for face in FW:
        pivotFace = [W[v] for v in face+(face[0],)]
        transform = submanifoldMapping(pivotFace)
        mappedVerts = (transform * (mat([p+[1.0] for p in pivotFace]).T)).T.tolist()
        facet = [point[:-2] for point in mappedVerts]
        pol = PolygonTessellator()
        vertices = [ vertex.Vertex( (x,y,0) ) for (x,y) in facet  ]
        verts = pol.tessellate(vertices)
        ps = [list(v.point) for v in verts]
        trias = [[ps[k],ps[k+1],ps[k+2],ps[k]] for k in range(0,len(ps),3)]
        mappedVerts = (transform.I * (mat([p+[1.0] for p in ps]).T)).T.tolist()
        points = [p[:-1] for p in mappedVerts]
        trias = [[points[k],points[k+1],points[k+2],points[k]] for k in range(0,len(points),3)]
        triangles += DISTR([AA(orientTriangle)(trias),[[0,1,2]]])
    return triangles

triangles = boundaryTriangulation(W,FW)
VIEW(EXPLODE(1.2,1.2,1.2)(CAT(AA(MKPOLS)(triangles))))
