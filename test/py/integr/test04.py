""" Integrals on 2D non-convex polyline """
from larlib import *

polyline = TRANS([[10,10,20,40,30,30,15,15],[10,20,30,20,10,15,15,10]])
model = polyline2lar([polyline])
VIEW(POLYLINE(polyline+[polyline[0]]))

V,FV,EV = model
tria = FV[0]
triangles = AA(C(AL)(0))(TRANS([tria[1:-1],tria[2:]]))
V = [v+[0.0] for v in V]
P = V,triangles

area = Surface(P,signed=True)
print "area =",area
