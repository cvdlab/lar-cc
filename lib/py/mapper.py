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


""" simplicial decomposition of the unit d-cube """
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
   points,verts = TRANS(vertDict.items())
   invertedindex = [None]*n
   V = []
   for k,value in enumerate(verts):
      V.append(eval(points[k]))
      for i in value:
         invertedindex[i]=k   
   CV = [[invertedindex[v] for v in cell] for cell in CV]
   # filter out degenerate cells
   CV = [list(set(cell)) for cell in CV if len(set(cell))==len(cell)]
   return V, CV

def larMap(coordFuncs):
   def larMap0(domain):
      V,CV = domain
      V = TRANS(CONS(coordFuncs)(V))  # plasm CONStruction
      return checkModel((V,CV))
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
      x = lambda V : [radius*COS(p[0]) for p in V]
      y = lambda V : [radius*SIN(p[0]) for p in V]
      return larMap([x,y])(domain)
   return larCircle0

def larDisk(radius=1.):
   def larDisk0(shape=[36,1]):
      domain = larIntervals(shape)([2*PI,radius])
      V,CV = domain
      x = lambda V : [p[1]*COS(p[0]) for p in V]
      y = lambda V : [p[1]*SIN(p[0]) for p in V]
      return larMap([x,y])(domain)
   return larDisk0

def larRing(params):
   r1,r2 = params
   def larRing0(shape=[36,1]):
      V,CV = larIntervals(shape)([2*PI,r2-r1])
      V = translatePoints(V,[0,r1])
      domain = V,CV
      x = lambda V : [p[1] * COS(p[0]) for p in V]
      y = lambda V : [p[1] * SIN(p[0]) for p in V]
      return larMap([x,y])(domain)
   return larRing0

def larSphere(radius=1):
   def larSphere0(shape=[18,36]):
      V,CV = larIntervals(shape)([PI,2*PI])
      V = translatePoints(V,[-PI/2,-PI])
      domain = V,CV
      x = lambda V : [radius*COS(p[0])*COS(p[1]) for p in V]
      y = lambda V : [radius*COS(p[0])*SIN(p[1]) for p in V]
      z = lambda V : [radius*SIN(p[0]) for p in V]
      return larMap([x,y,z])(domain)
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
      return larMap([x,y,z])(domain)
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
      return larMap([x,y,z])(domain)
   return larCrown0

def larBall(radius=1):
   def larBall0(shape=[18,36]):
      V,CV = checkModel(larSphere(radius)(shape))
      return V,[range(len(V))]
   return larBall0

def larRod(params):
   radius,height= params
   def larRod0(shape=[36,1]):
      V,CV = checkModel(larCylinder(params)(shape))
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
      return larMap([x,y,z])(domain)
   return larTorus0

def larPizza(params):
   r,R= params
   def larPizza0(shape=[24,36]):
      V,CV = checkModel(larCrown(params)(shape))
      return V,[range(len(V))]
   return larPizza0

def t(*args): 
   d = len(args)
   mat = scipy.identity(d+1)
   for k in range(d): 
      mat[k,d] = args[k]
   return mat.view(Mat)

def s(*args): 
   d = len(args)
   mat = scipy.identity(d+1)
   for k in range(d): 
      mat[k,k] = args[k]
   return mat.view(Mat)

def r(*args): 
   args = list(args)
   n = len(args)
   if n == 1: # rotation in 2D
      angle = args[0]; cos = COS(angle); sin = SIN(angle)
      mat = scipy.identity(3)
      mat[0,0] = cos;   mat[0,1] = -sin;
      mat[1,0] = sin;   mat[1,1] = cos;
   
   if n == 3: # rotation in 2D
      mat = scipy.identity(4)
      angle = VECTNORM(args); axis = UNITVECT(args)
      cos = COS(angle); sin = SIN(angle)
      if axis[1]==axis[2]==0.0:  # rotation about x
         mat[1,1] = cos;   mat[1,2] = -sin;
         mat[2,1] = sin;   mat[2,2] = cos;
      elif axis[0]==axis[2]==0.0:   # rotation about y
         mat[0,0] = cos;   mat[0,2] = sin;
         mat[2,0] = -sin;  mat[2,2] = cos;
      elif axis[0]==axis[1]==0.0:   # rotation about z
         mat[0,0] = cos;   mat[0,1] = -sin;
         mat[1,0] = sin;   mat[1,1] = cos;
      
      else:    # general 3D rotation (Rodrigues' rotation formula)   
         I = scipy.identity(3) ; u = axis
         Ux = scipy.array([
            [0,      -u[2],    u[1]],
            [u[2],      0,    -u[0]],
            [-u[1],   u[0],      0]])
         UU = scipy.array([
            [u[0]*u[0], u[0]*u[1],  u[0]*u[2]],
            [u[1]*u[0], u[1]*u[1],  u[1]*u[2]],
            [u[2]*u[0], u[2]*u[1],  u[2]*u[2]]])
         mat[:3,:3] = cos*I + sin*Ux + (1.0-cos)*UU
      
   
   return mat.view(Mat)

def larEmbed(k):
   def larEmbed0(model):
      V,CV = model
      if k>0:
         V = [v+[0.]*k for v in V] 
      elif k<0:
         V = [v[:-k] for v in V] 
      return V,CV
   return larEmbed0

def larApply(affineMatrix):
   def larApply0(model):
      if isinstance(model,Model):
         V = scipy.dot([v.tolist()+[1.0] for v in model.verts], affineMatrix.T).tolist()
         V = [v[:-1] for v in V]
         CV = copy(model.cells)
         d = copy(model.d)
         return Model((V,CV),d)
      elif isinstance(model,tuple):
         V,CV = model
         V = scipy.dot([v+[1.0] for v in V], affineMatrix.T).tolist()
         return [v[:-1] for v in V],CV
   return larApply0

""" Traversal of a scene multigraph """
def traverse(CTM, stack, o, scene=[]):
    for i in range(len(o)):
        if isinstance(o[i],Model): 
            scene += [larApply(CTM)(o[i])]
        elif isinstance(o[i],Mat): 
            CTM = scipy.dot(CTM, o[i])
        elif isinstance(o[i],Struct):
            stack.append(CTM) 
            traverse(CTM, stack, o[i], scene)
            CTM = stack.pop()
    return scene


def evalStruct(struct):
    dim = struct.n
    CTM, stack = scipy.identity(dim+1), []
    scene = traverse(CTM, stack, struct) 
    return scene

