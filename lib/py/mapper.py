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


