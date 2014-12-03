"""Module with automatic generation of simplified 3D buildings"""
import sys; sys.path.insert(0, 'lib/py/')
from architectural import *
""" transform svg primitives to basic lar format """
def rects2polylines(rooms): 
   return [[[x,y],[x+dx,y],[x+dx,y+dy],[x,y+dy]] for x,y,dx,dy in rooms]
   
def polyline2lar(polylines):
   index,defaultValue = -1,-1
   Vdict,FV = dict(),[]
   for k,polyline in enumerate(polylines):
      cell = []
      for vert in polyline:
         key = vcode(vert)
         if Vdict.get(key,defaultValue) == defaultValue:
            index += 1
            Vdict[key] = index
            cell += [index]
         else: 
            cell += [Vdict[key]]
      FV += [cell]
   items = TRANS(Vdict.items())
   V = TRANS(sorted(zip(items[1],AA(eval)(items[0]))))[1]
   #FV = AA(sorted)(FV)
   return V,FV

""" transform a lar model to a list of lar structures """
def lar2Structs(model):
   V,FV = model
   return [ Struct([[[V[v] for v in cell], [range(len(cell))]]]) for cell in FV]

""" transform an absolute lar model to a relative lar structure """
def absModel2relStruct(larPolylineModel):
   V,E = larPolylineModel
   Vnew = (array(V) - V[0]).tolist()
   return Struct([ t(*V[0]), (Vnew,E) ])

""" print a lar structure to a geoJson file """
def printStruct2GeoJson(struct):
   dim = checkStruct(struct.body)
   print "\n dim =",dim
   CTM, stack = scipy.identity(dim+1), []
   print "\n CTM, stack =",CTM, stack
   scene = printTraversal(CTM, stack, struct, [], 0) 
   return scene

""" Traverse a structure to print a geoJson file """
def printTraversal(CTM, stack, obj, scene=[], level=0):
   tabs = (4*level)*" "
   for i in range(len(obj)):
      if isinstance(obj[i],Model): 
         i,verts,cells = obj[i]
         if obj[i].__name__() == None:
            name = id(obj[i])
         else: 
            name = obj[i].__name__()
         print tabs, "i =",i
         print tabs, "name =",name
         print tabs, "verts =",AA(eval)(AA(vcode)(verts))
         print tabs, "cells =",cells
         scene += [larApply(CTM)(obj[i])]
      elif (isinstance(obj[i],tuple) or isinstance(obj[i],list)) and len(obj[i])==2:
         verts,cells = obj[i]
         name = id(obj[i])
         print tabs, "i =",i
         print tabs, "name =",name
         print tabs, "verts =",AA(eval)(AA(vcode)(verts))
         print tabs, "cells =",cells
         scene += [larApply(CTM)(obj[i])]
      elif isinstance(obj[i],Mat): 
         print tabs, "tVector =", obj[i].T[-1].tolist()
         CTM = scipy.dot(CTM, obj[i])
      elif isinstance(obj[i], Struct):
         if obj[i].__name__() == None:
            name = id(obj[i])
         else: 
            name = obj[i].__name__()
         print tabs, "i =",i
         print tabs, "name =",name
         stack.append(CTM) 
         level += 1
         printTraversal(CTM, stack, obj[i], scene, level)
         level -= 1
         CTM = stack.pop()
   return scene

