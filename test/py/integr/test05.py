""" Integrals on 2D non-convex polyline """
from larlib import *

V,FV,EV = openCourt11.body[0]
tria = FV[0]
triangles = AA(C(AL)(0))(TRANS([tria[1:-1],tria[2:]]))
V = [v+[0.0] for v in V]
P = V,triangles

area = Surface(P,signed=True)
print "area =",area

P = V, AA(REVERSE)(triangles)
area = Surface(P,signed=True)
print "area =",area
