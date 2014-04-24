"""Module with functions needed to interface LAR with pyplasm"""
import sys
sys.path.insert(0, 'lib/py/')
import simplexn

def importModule(moduleName):
   import sys
   sys.path.insert(0, 'lib/py/')
   import moduleName
   

import scipy as sp
from pyplasm import *
def CCOMB(vectors):
    return (sp.array(VECTSUM(vectors)) / float(len(vectors))).tolist()  

""" class definitions for LAR """
import scipy
class Mat(scipy.ndarray): pass
class Verts(scipy.ndarray): pass

class Model:
   """ A pair (geometry, topology) of the LAR package """
   def __init__(self,(verts,cells),dim):
      self.n = len(verts[0])
      self.d = dim
      self.verts = scipy.array(verts).view(Verts)
      self.cells = cells

def checkStruct(data):  
   def visit(o):
       if isinstance(o, (list,Struct)):
           for value in o:
               for subvalue in visit(value):
                   yield subvalue
       elif isinstance(o, Model) or isinstance(o, Mat): 
         yield o     
   flatten = list(visit(data))
   dims = [ o.n for o in flatten if isinstance(o, Model) ]
   if EQ(dims): return dims[0] 
   else:return None

class Struct:
    """ The assembly type of the LAR package """
    def __init__(self,data):
        self.n = checkStruct(data)
        self.body = data
    def __iter__(self):
        return iter(self.body)
    def __len__(self):
        return len(list(self.body))
    def __getitem__(self,i):
        return list(self.body)[i]

def MKPOLS (model):
    if isinstance(model,Model):
        V, FV = model.verts.tolist(), model.cells
    elif isinstance(model,tuple) or isinstance(model,list):
        V, FV = model
    pols = [MKPOL([[V[v] for v in f],[range(1,len(f)+1)], None]) for f in FV]
    return pols  

def EXPLODE (sx,sy,sz):
    def explode0 (scene):
        centers = [CCOMB(S1(UKPOL(obj))) for obj in scene]
        scalings = len(centers) * [S([1,2,3])([sx,sy,sz])]
        scaledCenters = [UK(APPLY(pair)) for pair in
                         zip(scalings, [MK(p) for p in centers])]
        translVectors = [ VECTDIFF((p,q)) for (p,q) in zip(scaledCenters, centers) ]
        translations = [ T([1,2,3])(v) for v in translVectors ]
        return STRUCT([ t(obj) for (t,obj) in zip(translations,scene) ])
    return explode0  

