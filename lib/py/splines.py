""" Mapping functions and primitive objects """
from pyplasm import *
from scipy import *
import os,sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from lar2psm import *
from simplexn import *
from larcc import *
from largrid import *
from mapper import *

""" Multidimensional transfinite Bezier """
def larBezier(U,d=3):
   def BEZIER0(controldata_fn):
      N = len(controldata_fn)-1
      def map_fn(point):
         t = U(point)
         controldata = [fun(point) if callable(fun) else fun 
            for fun in controldata_fn]
         out = [0.0 for i in range(len(controldata[0]))]    
         for I in range(N+1):
            weight = CHOOSE([N,I])*math.pow(1-t,N-I)*math.pow(t,I)
            for K in range(len(out)):  out[K] += weight*(controldata[I][K])
         return out
      return (COMP([AA(COMP),DISTR]))([AA(SEL)(range(d)), map_fn])
   return BEZIER0

def larBezierCurve(controlpoints):
   dim = len(controlpoints[0])
   return larBezier(S1,dim)(controlpoints)

""" Transfinite Coons patches """
def larCoonsPatch (args):
   su0_fn , su1_fn , s0v_fn , s1v_fn = args
   def map_fn(point):
      u,v=point
      su0 = su0_fn(point) if callable(su0_fn) else su0_fn
      su1 = su1_fn(point) if callable(su1_fn) else su1_fn
      s0v = s0v_fn(point) if callable(s0v_fn) else s0v_fn
      s1v = s1v_fn(point) if callable(s1v_fn) else s1v_fn
      ret=[0.0 for i in range(len(su0))]  
      for K in range(len(ret)):
         ret[K] = ((1-u)*s0v[K] + u*s1v[K]+(1-v)*su0[K] + v*su1[K] + 
         (1-u)*(1-v)*s0v[K] + (1-u)*v*s0v[K] + u*(1-v)*s1v[K] + u*v*s1v[K])
      return ret
   return (COMP([AA(COMP),DISTR]))([[S1,S2,S3], map_fn])

