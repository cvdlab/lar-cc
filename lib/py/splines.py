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

""" Tensor product surface patch """
def larTensorProdSurface (args):
   ubasis , vbasis = args
   def TENSORPRODSURFACE0 (controlpoints_fn):
      def map_fn(point):
         u,v=point
         U=[f([u]) for f in ubasis]
         V=[f([v]) for f in vbasis]
         controlpoints=[f(point) if callable(f) else f 
            for f in controlpoints_fn]
         target_dim = len(controlpoints[0][0])
         ret=[0 for x in range(target_dim)]
         for i in range(len(ubasis)):
            for j in range(len(vbasis)):
               for M in range(len(ret)):
                  for M in range(target_dim): 
                     ret[M] += U[i]*V[j] * controlpoints[i][j][M]
         return ret
      return map_fn
   return TENSORPRODSURFACE0

""" Bilinear tensor product surface patch """
def larBilinearSurface(controlpoints):
   basis = larBernsteinBasis(S1)(1)
   return larTensorProdSurface([basis,basis])(controlpoints)

""" Biquadratic tensor product surface patch """
def larBiquadraticSurface(controlpoints):
   basis1 = larBernsteinBasis(S1)(2)
   basis2 = larBernsteinBasis(S1)(2)
   return larTensorProdSurface([basis1,basis2])(controlpoints)

""" Bicubic tensor product surface patch """
def larBicubicSurface(controlpoints):
   basis1 = larBernsteinBasis(S1)(3)
   basis2 = larBernsteinBasis(S1)(3)
   return larTensorProdSurface([basis1,basis2])(controlpoints)

""" Toolbox of tensor operations """
def larBernsteinBasis (U):
   def BERNSTEIN0 (N):
      def BERNSTEIN1 (I):
         def map_fn(point):
            t = U(point)
            out = CHOOSE([N,I])*math.pow(1-t,N-I)*math.pow(t,I)
            return out
         return map_fn
      return [BERNSTEIN1(I) for I in range(0,N+1)]
   return BERNSTEIN0

""" Multidimensional transfinite Bezier """
def larBezier(U):
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
      return map_fn
   return BEZIER0

def larBezierCurve(controlpoints):
   return larBezier(S1)(controlpoints)

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
   return map_fn

""" Domain decomposition for 1D bspline maps """
def larDom(knots,tics=32): 
   domain = knots[-1]-knots[0]
   return larIntervals([tics*domain])([domain])

""" Alias for the pyplasm definition (too long :o) """
NURBS = RATIONALBSPLINE

