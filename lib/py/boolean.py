""" Module with Boolean operators using chains and CSR matrices """
from pyplasm import *
from scipy import *
import os,sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from lar2psm import *
from simplexn import *
from larcc import *
from largrid import *
from myfont import *
from mapper import *

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

""" High level Boolean Application Programming Interface """
def larUnion(lar1,lar2): lar = boolOps(lar1,lar2); pass
def larIntersection(lar1,lar2): lar = boolOps(lar1,lar2); pass
def larDifference(lar1,lar2): lar = boolOps(lar1,lar2); pass
def larXor(lar1,lar2): lar = boolOps(lar1,lar2); pass

""" First step of Boolean Algorithm """
from collections import defaultdict, OrderedDict

""" TODO: change defaultdict to OrderedDefaultdict """

class OrderedDefaultdict(collections.OrderedDict):
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or callable(args[0])):
                raise TypeError('first argument must be callable or None')
            self.default_factory = args[0]
            args = args[1:]
        super(OrderedDefaultdict, self).__init__(*args, **kwargs)

    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = default = self.default_factory()
        return default

    def __reduce__(self):  # optional, for pickle support
        args = (self.default_factory,) if self.default_factory else tuple()
        return self.__class__, args, None, None, self.iteritems()


def vertexSieve(model1, model2):
   from lar2psm import larModelBreak
   V1,CV1 = larModelBreak(model1) 
   V2,CV2 = larModelBreak(model2)
   n = len(V1); m = len(V2)
   def shift(CV, n): 
      return [[v+n for v in cell] for cell in CV]
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


def boolOps(lar1,lar2,cell='simplex', facets1=None,facets2=None):
   (V1,CV1),(V2,CV2) = lar1,lar2
   n1,n2 = len(V1),len(V2)
   V, CV1, CV2, n12 = vertexSieve(lar1, lar2)
   CV = Delaunay(array(V)).vertices
   BV1, BV2 = boundaryVertices( V, CV1,CV2, cell, facets1,facets2 )
   print "\n BV1 =",BV1
   print "\n BV2 =",BV2
   """ Delaunay triangulation of boundary vertices """
   B = [V[v] for v in BV1+BV2]
   CV = Delaunay(array(B)).vertices
   VIEW(STRUCT([
      EXPLODE(1.2,1.2,1.2)(MKPOLS((B,CV))),
      COLOR(CYAN)(SKEL_1(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,CV1))))),
      COLOR(MAGENTA)(SKEL_1(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,CV2)))))
   ]))
   """ back-indicize globally (original common vertices) """
   BV = [v for k,v in enumerate(BV1+BV2)]
   CV = [[BV[v] for v in cell] for cell in CV]
   
   """ Common simplices extraction """
   CV_un, CV_int = splitDelaunayComplex(CV,n1,n2,n12)
   print "\n CV_un =",CV_un
   print "\n CV_int =",CV_int
   VIEW(COLOR(YELLOW)(EXPLODE(1.2,1.2,1)(MKPOLS((V,CV_int)))))
   
   return V,n1,n2,n12, BV1, BV2

""" Second stage of Boolean operations """
def boundaryVertices( V, CV1,CV2, cell='simplex', facets1=None,facets2=None ):
   if cell=='simplex': 
      FV1 = larSimplexFacets(CV1)
      FV2 = larSimplexFacets(CV2)
   elif cell=='cuboid': 
      FV1 = facets1
      print "\n FV1 =",FV1
      FV2 = facets2
      print "\n FV2 =",FV2
   BF1 = boundaryCells(CV1,FV1)
   print "\n BF1 =",BF1
   BF2 = boundaryCells(CV2,FV2)
   BV1 = list(set(CAT([ FV1[f] for f in BF1 ])))
   BV2 = list(set(CAT([ FV2[f] for f in BF2 ])))
   VIEW(STRUCT([ 
      COLOR(GREEN)(STRUCT(AA(MK)([V[v] for v in BV1]))), 
      COLOR(YELLOW)(STRUCT(AA(MK)([V[v] for v in BV2]))) ]))
   return BV1, BV2

def splitDelaunayComplex(CV,n1,n2,n12):
   def test(cell):
      return any([v<n1 for v in cell]) and any([v>=(n1-n12) for v in cell])
   cells_intersection, cells_union = [],[]
   for cell in CV: 
      if test(cell): cells_intersection.append(cell)
      else: cells_union.append(cell)
   return cells_union,cells_intersection

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


def cellNames(model,cells, color=BLACK):
   V,CV= model
   print "\n CV =",CV
   print "\n cells =",cells
   texts = []
   for k,cell in enumerate(cells):
      centroid = CCOMB([V[v] for v in cell])
      print "centroid =",centroid
      d = len(centroid)
      texts += [ T(range(1,d+1))(centroid)(S(range(1,d+1))([0.02 
                  for h in range(d)])(TEXTWITHATTRIBUTES()(str(k)))) ]
   return AA(COLOR(color))(texts)

