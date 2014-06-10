""" Boundary orientation of a random 2D triangulation """
import sys;sys.path.insert(0, 'lib/py/')
from scipy import linalg
from larcc import *
from random import random

""" Vertices V generated as random point in the unit circle """
verts = []
npoints = 200
for k in range(npoints):
   t = 2*pi*random()
   u = random()+random()
   if u > 1: r = 2-u 
   else: r = u
   verts += [[r*cos(t), r*sin(t)]]
VIEW(STRUCT(AA(MK)(verts)))

""" Delaunay triangulation of the whole set V of points """
triangles = Delaunay(verts)
def area(cell): return linalg.det([verts[v]+[1] for v in cell])/2
cells = [ cell for cell in triangles.vertices.tolist() if area(cell)>PI/(3*npoints)]
V, FV = AA(list)(verts), cells

""" Fraction of triangles randomly discarded """
fraction = 0.7
cellSpan = len(FV)
remove = [int(random()*cellSpan) for k in range(int(cellSpan*fraction)) ]
FV = [FV[k] for k in range(cellSpan) if not k in remove]

""" Coherently orient the input LAR model (V,FV) """
def positiveOrientation(model):
   V,simplices = model
   out = []
   for simplex in simplices:
      theMat = [V[v]+[1] for v in simplex]
      if sign(linalg.det(theMat)) > 0:  out += [simplex]
      else: out += [REVERSE(simplex)]
   return V,out

V,FV = positiveOrientation((V,FV))

""" Compute the 1-cell and 0-cell bases EV and VV """
EV = larSimplexFacets(FV)
VV = AA(LIST)(range(len(V)))
VIEW(mkSignedEdges((V,EV)))

""" Signed 2-boundary matrix  and signed boundary 1-chain """
orientedBoundary = signedCellularBoundaryCells(V,[VV,EV,FV])

""" Display the boundary 1-chain """
VIEW(STRUCT(MKPOLS((V,FV))))
VIEW(STRUCT(
   MKPOLS((V,FV)) +
   [COLOR(RED)(mkSignedEdges((V,orientedBoundary)))]  ))

