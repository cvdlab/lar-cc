""" Mapping functions and primitive objects """
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

import boolean2
from boolean2 import *


def larDomain(shape):
   V,CV = larSimplexGrid(shape)
   V = scalePoints(V, [1./d for d in shape])
   return V,CV

def larIntervals(shape):
   def larIntervals0(size):
      V,CV = larDomain(shape)
      V = scalePoints(V, [scaleFactor for scaleFactor in size])
      return V,CV
   return larIntervals0

def checkModel(model):
   V,CV = model; n = len(V)
   vertDict = defaultdict(list)
   for k,v in enumerate(V): vertDict[vcode(v)].append(k) 
   verts = (vertDict.values())
   invertedindex = [None]*n
   for k,value in enumerate(verts):
      print "\n len(value) =",len(value)
      for i in value:
         invertedindex[i]=value[0]  
   CV = [[invertedindex[v] for v in cell] for cell in CV]
   # filter out degenerate cells
   CV = [list(set(cell)) for cell in CV if len(set(cell))==len(cell)]
   return V, CV

def larMap(coordFuncs):
   def larMap0(domain):
      V,CV = domain
      V = TRANS(CONS(coordFuncs)(V))
      return V,CV
   return larMap0

if __name__=="__main__":
   V,EV = larDomain([5])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,EV))))
   V,EV = larIntervals([24])([2*PI])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,EV))))
      
   V,FV = larDomain([5,3])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))
   V,FV = larIntervals([36,3])([2*PI,1.])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))
      
   V,CV = larDomain([5,3,1])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))
   V,CV = larIntervals([36,2,3])([2*PI,1.,1.])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))

def larCircle(radius=1.):
   def larCircle0(shape=36):
      domain = larIntervals([shape])([2*PI])
      V,CV = domain
      x = lambda coords : [radius*COS(p[0]) for p in V]
      y = lambda coords : [radius*SIN(p[0]) for p in V]
      mapping = [x,y]
      return larMap(mapping)(domain)
   return larCircle0

def larDisk(radius=1.):
   def larDisk0(shape=[36,1]):
      domain = larIntervals(shape)([2*PI,radius])
      V,CV = domain
      x = lambda V : [p[1]*COS(p[0]) for p in V]
      y = lambda V : [p[1]*SIN(p[0]) for p in V]
      mapping = [x,y]
      return larMap(mapping)(domain)
   return larDisk0


def larRing(params):
   r1,r2 = params
   def larDisk0(shape=[36,1]):
      V,CV = larIntervals(shape)([2*PI,r2-r1])
      V = translatePoints(V,[0,r1])
      domain = V,CV
      VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,CV))))
      x = lambda V : [p[1] * COS(p[0]) for p in V]
      y = lambda V : [p[1] * SIN(p[0]) for p in V]
      mapping = [x,y]
      return larMap(mapping)(domain)
   return larDisk0


def larSphere(radius=1):
   def larSphere0(shape=[18,36]):
      V,CV = larIntervals(shape)([PI,2*PI])
      V = translatePoints(V,[-PI/2,-PI])
      domain = V,CV
      x = lambda V : [radius*COS(p[0])*SIN(p[1]) for p in V]
      y = lambda V : [radius*COS(p[0])*COS(p[1]) for p in V]
      z = lambda V : [radius*SIN(p[0]) for p in V]
      mapping = [x,y,z]
      return larMap(mapping)(domain)
   return larSphere0

from scipy.linalg import det
"""
def makeOriented(model):
   V,CV = model
   out = []
   for cell in CV: 
      mat = scipy.array([V[v]+[1] for v in cell]+[[0,0,0,1]])
      if det(mat) < 0.0:
         out.append(cell)
      else:
         out.append([cell[1]]+[cell[0]]+cell[2:])
      print "\n det(mat) =",det(mat)
   return V,out
"""
def larCylinder(params):
   radius,height= params
   def larCylinder0(shape=[36,1]):
      domain = larIntervals(shape)([2*PI,1])
      V,CV = domain
      x = lambda V : [radius*COS(p[0]) for p in V]
      y = lambda V : [radius*SIN(p[0]) for p in V]
      z = lambda V : [height*p[1] for p in V]
      mapping = [x,y,z]
      model = larMap(mapping)(domain)
      # model = makeOriented(model)
      return model
   return larCylinder0

def larToroidal(params):
   r,R = params
   def larToroidal0(shape=[24,36]):
      domain = larIntervals(shape)([2*PI,2*PI])
      V,CV = domain
      x = lambda V : [(R + r*COS(p[0])) * COS(p[1]) for p in V]
      y = lambda V : [(R + r*COS(p[0])) * SIN(p[1]) for p in V]
      z = lambda V : [-r * SIN(p[0]) for p in V]
      mapping = [x,y,z]
      return larMap(mapping)(domain)
   return larToroidal0

def larCrown(params):
   r,R = params
   def larCrown0(shape=[24,36]):
      V,CV = larIntervals(shape)([PI,2*PI])
      V = translatePoints(V,[-PI/2,0])
      domain = V,CV
      x = lambda V : [(R + r*COS(p[0])) * COS(p[1]) for p in V]
      y = lambda V : [(R + r*COS(p[0])) * SIN(p[1]) for p in V]
      z = lambda V : [-r * SIN(p[0]) for p in V]
      mapping = [x,y,z]
      return larMap(mapping)(domain)
   return larCrown0

def larBall(radius=1):
   def larBall0(shape=[18,36]):
      V,CV = checkModel(larSphere(radius)(shape))
      VIEW(STRUCT(MKPOLS((V,CV))))
      return V,[range(len(V))]
   return larBall0

def larRod(params):
   radius,height= params
   def larRod0(shape=[36,1]):
      V,CV = checkModel(larCylinder(params)(shape))
      VIEW(STRUCT(MKPOLS((V,CV))))
      return V,[range(len(V))]
   return larRod0

def larTorus(params):
   r,R = params
   def larTorus0(shape=[24,36,1]):
      domain = larIntervals(shape)([2*PI,2*PI,r])
      V,CV = domain
      x = lambda V : [(R + p[2]*COS(p[0])) * COS(p[1]) for p in V]
      y = lambda V : [(R + p[2]*COS(p[0])) * SIN(p[1]) for p in V]
      z = lambda V : [-p[2] * SIN(p[0]) for p in V]
      mapping = [x,y,z]
      return larMap(mapping)(domain)
   return larTorus0

def larPizza(params):
   r,R= params
   def larPizza0(shape=[24,36]):
      V,CV = checkModel(larCrown(params)(shape))
      VIEW(STRUCT(MKPOLS((V,CV))))
      return V,[range(len(V))]
   return larPizza0

