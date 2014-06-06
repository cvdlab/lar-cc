"""Module with functions needed to interface LAR with pyplasm"""
def importModule(moduleName):
   import sys; sys.path.insert(0, 'lib/py/')
   import moduleName
   

import scipy as sp
from pyplasm import *
def CCOMB(vectors):
    return (sp.array(VECTSUM(vectors)) / float(len(vectors))).tolist()  

import simplexn
from simplexn import *
""" class definitions for LAR """
import scipy
class Mat(scipy.ndarray): pass
class Verts(scipy.ndarray): pass

class Model:
   """ A pair (geometry, topology) of the LAR package """
   def __init__(self,(verts,cells)):
      self.n = len(verts[0])
      # self.verts = scipy.array(verts).view(Verts)
      self.verts = verts
      self.cells = cells

class Struct:
    """ The assembly type of the LAR package """
    def __init__(self,data):
        self.body = data
    def __iter__(self):
        return iter(self.body)
    def __len__(self):
        return len(list(self.body))
    def __getitem__(self,i):
        return list(self.body)[i]

def larModelBreak(model):
    if isinstance(model,Model):
        # V, FV = model.verts.tolist(), model.cells
        V, FV = model.verts, model.cells
    elif isinstance(model,tuple) or isinstance(model,list):
        V, FV = model
    return V,FV

def MKPOLS (model):
   V,FV = larModelBreak(model)
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

