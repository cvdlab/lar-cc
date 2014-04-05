""" Module with Boolean operators using chains and CSR matrices """
from pyplasm import *
from scipy import *
import os,sys

""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
import lar2psm
from lar2psm import *

import simplexn
from simplexn import *

import larcc
from larcc import *

import largrid
from largrid import *

import myfont
from myfont import *


""" TODO: use package Decimal (http://docs.python.org/2/library/decimal.html) """
ROUND_ZERO = 1E-07
def round_or_zero (x,prec=7):
   """
   Decision procedure to approximate a small number to zero.
   Return either the input number or zero.
   """
   def myround(x):
      return eval(('%.'+str(prec)+'f') % round(x,prec))
   xx = myround(x)
   if abs(xx) < ROUND_ZERO: return 0.0
   else: return xx

def prepKey (args): return "["+", ".join(args)+"]"

def fixedPrec(value):
   if abs(value - int(value))<ROUND_ZERO: value = int(value)
   out = ('%0.7f'% value).rstrip('0')
   if out == '-0.': out = '0.'
   return out
   
def vcode (vect): 
   """
   To generate a string representation of a number array.
   Used to generate the vertex keys in PointSet dictionary, and other similar operations.
   """
   return prepKey(AA(fixedPrec)(vect))

def translatePoints (points, tvect):
   return [VECTSUM([p,tvect]) for p in points]

def rotatePoints (points, angle):
   return [[COS(x),-SIN(y)] for x,y in points]

def scalePoints (points, svect):
   return [AA(PROD)(TRANS([p,svect])) for p in points]

def randomPointsInUnitCircle(n=200,d=2, r=1):
   points = random.random((n,d)) * ([2*math.pi]+[1]*(d-1))
   return [[SQRT(p[1])*COS(p[0]),SQRT(p[1])*SIN(p[0])] for p in points]
   ## TODO: correct for $d$-sphere

if __name__=="__main__":
   VIEW(STRUCT(AA(MK)(randomPointsInUnitCircle()))) 

def randomPointsInUnitCuboid(n=200,d=2):
   return random.random((n,d)).tolist()

if __name__=="__main__":
   VIEW(STRUCT(AA(MK)(randomPointsInUnitCuboid()))) 

from scipy.spatial import Delaunay
def randomTriangulation(n=200,d=2,out='disk'):
   if out == 'disk':
      V = randomPointsInUnitCircle(n,d)
   elif out == 'cuboid':
      V = randomPointsInUnitCuboid(n,d)
   CV = Delaunay(array(V)).vertices
   model = V,CV
   return model

if __name__=="__main__":
   from lar2psm import *
   VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(model)))

""" High level Boolean Application Programming Interface """

def boolOps(lar1,lar2):
   (V1,CV1),(V2,CV2) = lar1,lar2
   n1,n2 = len(V1),len(V2)
   V, CV1, CV2, n12 = vertexSieve(lar1, lar2)
   CV = Delaunay(array(V)).vertices
   setCV = set([tuple(sorted(cell)) for cell in CV])
   # Second stage of Boolean algorithm
   B1,B2 = boundaryVertices( V, CV1, CV2 )
   # Extraction of coboundary of boundary chains
   BSupCells1 = boundarySuperCells( V, CV1 )
   BSupCells2 = boundarySuperCells( V, CV2 )
   
   VIEW(STRUCT([ 
      COLOR(GREEN)(STRUCT(MKPOLS(((V,[CV1[c] for c in BSupCells1]))))), 
      COLOR(MAGENTA)(STRUCT(MKPOLS(((V,[CV2[c] for c in BSupCells2]))))) 
   ])) 
   
   def invariantSuperCells(V,BSupCells1,CV1,setCV):
      out = []
      for cell in BSupCells1:
         if tuple(CV1[cell]) in setCV: out += [cell]
      return out
   
   # Invariant subsets of supercells
   invSupCells1 = invariantSuperCells(V,BSupCells1,CV1,setCV)
   invSupCells2 = invariantSuperCells(V,BSupCells2,CV2,setCV)
   
   divisor1 = set(BSupCells1).difference(invSupCells1)
   divisor2 = set(BSupCells2).difference(invSupCells2)
   VIEW(STRUCT([ 
      COLOR(GREEN)(STRUCT(MKPOLS(((V,[CV1[c] for c in divisor1]))))), 
      COLOR(MAGENTA)(STRUCT(MKPOLS(((V,[CV2[c] for c in divisor2]))))) 
   ]))
   
   """ Minimal covering chains of divisors """
   def minimalCovers(V,divisors,CV1,CV,setCV):
      covers = [selectIncidentChain( V, CV )(v) for cell in divisors for v in CV1[cell] ]
      print "\n covers =",covers
      return covers
   
   
   # Covers of divisors computation
   covers1 = minimalCovers(V,divisor1,CV1,CV,setCV)
   covers2 = minimalCovers(V,divisor2,CV2,CV,setCV)   
   VIEW(STRUCT([ 
      EXPLODE(1.2,1.2,1)(MKPOLS(((V,CV)))), 
      COLOR(GREEN)(EXPLODE(1.2,1.2,1)(MKPOLS(((V,[CV[c] for c in CAT(covers1)]))))), 
      COLOR(MAGENTA)(EXPLODE(1.2,1.2,1)(MKPOLS(((V,[CV[c] for c in CAT(covers2)])))))
       ]))
   
   return V,n1,n2,n12, B1,B2

def union(lar1,lar2):
   lar = boolOps(lar1,lar2)
def intersection(lar1,lar2):
   lar = boolOps(lar1,lar2)
def difference(lar1,lar2):
   lar = boolOps(lar1,lar2)
def xor(lar1,lar2):
   lar = boolOps(lar1,lar2)

""" First step of Boolen Algorithm """
from collections import defaultdict, OrderedDict

def vertexSieve(model1, model2):
   V1,CV1 = model1; V2,CV2 = model2
   n = len(V1); m = len(V2)
   def shift(CV, n): 
      return [[v+n for v in cell]for cell in CV]
   CV2 = shift(CV2,n)


   vdict1 = defaultdict(list)
   for k,v in enumerate(V1): vdict1[vcode(v)].append(k) 
   vdict2 = defaultdict(list)
   for k,v in enumerate(V2): vdict2[vcode(v)].append(k+n) 
   
   vertdict = defaultdict(list)
   for point in vdict1.keys(): vertdict[point] += vdict1[point]
   for point in vdict2.keys(): vertdict[point] += vdict2[point]


   case1, case12, case2 = [],[],[]
   for item in vertdict.items():
      key,val = item
      if len(val)==2:  case12 += [item]
      elif val[0] < n: case1 += [item]
      else: case2 += [item]
   n1 = len(case1); n2 = len(case12); n3 = len(case2)


   invertedindex = list(0 for k in range(n+m))
   for k,item in enumerate(case1):
      invertedindex[item[1][0]] = k
   for k,item in enumerate(case12):
      invertedindex[item[1][0]] = k+n1
      invertedindex[item[1][1]] = k+n1
   for k,item in enumerate(case2):
      invertedindex[item[1][0]] = k+n1+n2


   V = [eval(p[0]) for p in case1] + [eval(p[0]) for p in case12] + [eval(
            p[0]) for p in case2]
   CV1 = [sorted([invertedindex[v] for v in cell]) for cell in CV1]
   CV2 = [sorted([invertedindex[v] for v in cell]) for cell in CV2]
   return V, CV1, CV2, len(case12)



""" Second stage of Boolean operations """
def boundaryVertices( V, CV1,CV2 ):
   FV1 = larSimplexFacets(CV1)
   FV2 = larSimplexFacets(CV2)
   BF1 = boundaryCells(CV1,FV1)
   BF2 = boundaryCells(CV2,FV2)
   BV1 = list(set([ v for f in BF1 for v in FV1[f] ]))
   BV2 = list(set([ v for f in BF2 for v in FV2[f] ]))
   VIEW(STRUCT([ 
      COLOR(GREEN)(STRUCT(AA(MK)([V[v] for v in BV1]))), 
      COLOR(MAGENTA)(STRUCT(AA(MK)([V[v] for v in BV2]))) ]))
   return BV1, BV2

""" Select the $d$-chain incident on a $0$-chain """
def selectIncidentChain( V, CV ):
   def selectIncidentChain0( v ):
      return [k for k in range(len(CV)) if v in CV[k] ]
   return selectIncidentChain0

def coboundaryOfBoundaryCells(cells,facets):
    csrBoundaryMat = boundary(cells,facets)
    csrChain = totalChain(cells)
    csrBoundaryChain = matrixProduct(csrBoundaryMat, csrChain)
    for k,value in enumerate(csrBoundaryChain.data):
        if value % 2 == 0: csrBoundaryChain.data[k] = 0
    boundaryCells = [k for k,val in enumerate(csrBoundaryChain.data.tolist())
                               if val == 1]
    csrCoboundaryBoundaryChain = matrixProduct(csrBoundaryMat.T, csrBoundaryChain)
    coboundaryBoundaryCells = [k for k in csrCoboundaryBoundaryChain.indices.tolist()]
    return coboundaryBoundaryCells

""" Second stage of Boolean operations """
def boundarySuperCells( V, CV ):
   FV = larSimplexFacets(CV)
   BSupCells = coboundaryOfBoundaryCells(CV,FV)
   return BSupCells

def cellNames(model,cells, color=BLACK):
   V,CV= model
   print "\n CV =",CV
   print "\n cells =",cells
   texts = []
   for k,cell in enumerate(cells):
      centroid = CCOMB([V[v] for v in cell])
      print "centroid =",centroid
      d = len(centroid)
      texts += [T(range(1,d+1))(centroid)(S(range(1,d+1))([0.005 for k in range(d)])(TEXT(str(k))))]
   return AA(COLOR(color))(texts)

