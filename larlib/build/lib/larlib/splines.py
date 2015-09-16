""" Mapping functions and primitive objects """
from larlib import *

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
   return larIntervals([tics*int(domain)])([domain])

""" Sampling of a set of B-splines of given degree, knots and controls """
def BSPLINEBASIS(degree):
   def BSPLINE0(knots):
      def BSPLINE1(ncontrols):
         n = ncontrols-1
         m=len(knots)-1
         k=degree+1
         T=knots
         tmin,tmax=T[k-1],T[n+1]       
         if len(knots)!=(n+k+1):
            raise Exception("Invalid point/knots/degree for bspline!")        

         # de Boor coefficients
         def N(i,k,t):           
            # Ni1(t)
            if k==1:
               if(t>=T[i] and t<T[i+1]) or(t==tmax and t>=T[i] and t<=T[i+1]):
                  # i use strict inclusion for the max value
                  return 1
               else:
                  return 0          
            # Nik(t)
            ret=0
            
            num1,div1= t-T[i], T[i+k-1]-T[i]
            if div1!=0: ret+=(num1/div1) * N(i,k-1,t)          
            num2,div2=T[i+k]-t, T[i+k]-T[i+1]
            if div2!=0:  ret+=(num2/div2) * N(i+1,k-1,t)
            
            return ret
         
         # map function
         def map_fn(point):
            t=point[0]
            return [N(i,k,t) for i in range(n+1)]
                     
         return map_fn
      return BSPLINE1
   return BSPLINE0

""" Drawing the graph of a set of B-splines """
if __name__=="__main__":

   knots = [0,0,0,1,1,2,2,3,3,4,4,4]
   ncontrols = 9
   degree = 2
   obj = larMap(BSPLINEBASIS(degree)(knots)(ncontrols))(larDom(knots))
   
   funs = TRANS(obj[0])
   var = AA(CAT)(larDom(knots)[0])
   cells = larDom(knots)[1]
   
   graphs =  [[TRANS([var,fun]),cells] for fun in funs]
   graph = STRUCT(CAT(AA(MKPOLS)(graphs)))
   VIEW(graph)
   VIEW(STRUCT(MKPOLS(graphs[0]) + MKPOLS(graphs[-1])))


def TBSPLINE(U):
   def TBSPLINE0(degree):
      def TBSPLINE1(knots):
         def TBSPLINE2(points_fn):
   
            n=len(points_fn)-1
            m=len(knots)-1
            k=degree+1
            T=knots
            tmin,tmax=T[k-1],T[n+1]
   
            # see http://www.na.iac.cnr.it/~bdv/cagd/spline/B-spline/bspline-curve.html
            if len(knots)!=(n+k+1):
               raise Exception("Invalid point/knots/degree for bspline!")
   
            # de boord coefficients
            def N(i,k,t):
   
               # Ni1(t)
               if k==1: 
                  if(t>=T[i] and t<T[i+1]) or (t==tmax and t>=T[i] and t<=T[i+1]): 
                     # i use strict inclusion for the max value
                     return 1
                  else:
                     return 0
   
               # Nik(t)
               ret=0
   
               num1,div1= t-T[i], T[i+k-1]-T[i]  
               if div1!=0: ret+=(num1/div1) * N(i,k-1,t)
               # elif num1!=0: ret+=N(i,k-1,t)
   
               num2,div2=T[i+k]-t, T[i+k]-T[i+1]
               if div2!=0:  ret+=(num2/div2) * N(i+1,k-1,t)
               # elif num2!=0: ret+=N(i,k-1,t)
   
               return ret
   
            # map function
            def map_fn(point):
               t=U(point)
   
               # if control points are functions
               points=[f(point) if callable(f) else f for f in points_fn]
   
               target_dim=len(points[0])
               ret=[0 for i in range(target_dim)];
               for i in range(n+1):
                  coeff=N(i,k,t) 
                  for M in range(target_dim):
                     ret[M]+=points[i][M]*coeff
               return ret
   
            return map_fn
   
         return TBSPLINE2
      return TBSPLINE1
   return TBSPLINE0

""" Transfinite NURBS """
def TRATIONALBSPLINE(U):
   def TRATIONALBSPLINE0(degree):
      def TRATIONALBSPLINE1(knots):
         def TRATIONALBSPLINE2(points):
            bspline=TBSPLINE(U)(degree)(knots)(points)
            def map_fn(point):         
               ret=bspline(point)
               last=ret[-1]
               if last!=0: ret=[value/last for value in ret]
               ret=ret[:-1]
               return ret
            return map_fn
         return TRATIONALBSPLINE2
      return TRATIONALBSPLINE1
   return TRATIONALBSPLINE0

""" Alias for the pyplasm definition (too long :o) """
NURBS = RATIONALBSPLINE     # in pyplasm
TNURBS = TRATIONALBSPLINE   # in lar-cc (only)

