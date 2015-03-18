"""Module with automatic generation of simplified 3D buildings"""
import sys; sys.path.insert(0, 'lib/py/')
from architectural import *
class Struct2(Struct):
   def flatten(self): 
      structs = copy(self.body)
      structList = lar2Structs(CAT(structs.body))
      return [ absModel2relStruct(struct[0].body) for struct in structList ]

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
   EV = face2edge(FV)
   return V,FV,EV

""" transform a lar model to a list of lar structures """
def lar2Structs(model):
   V,FV,_ = model
   return [ Struct([[[V[v] for v in cell], [range(len(cell))]]]) for cell in FV]

""" transform an absolute lar model to a relative lar structure """
def absModel2relStruct(larModel):
   V,E = larModel
   Vnew = (array(V) - V[0]).tolist()
   return Struct([ t(*V[0]), (Vnew,E) ])

""" print a lar structure to a geoJson file """
import yaml
import json

def printStruct2GeoJson(path,struct):
    if struct.__name__() == None:
        filename = str(id(struct))
    else: 
        filename = struct.__name__()
    theFile = open(path+filename+".yml", "w")
    print >> theFile, "---"
    print "filename =", path+filename+".yml"
    dim = checkStruct(struct.body)
    CTM, stack = scipy.identity(dim+1), []
    fathers = [filename]
    scene,fathers = printTraversal(theFile, CTM, stack, struct, [], 0, fathers,filename) 
    theFile.close()
    """ exporting to JSON via YML a lar structure """
    file_handle = open(path+filename+".yml")
    my_dictionary = yaml.safe_load(file_handle)
    file_handle.close()
    with open(path+filename+".json", 'w') as outfile:
         json.dump(my_dictionary, outfile, sort_keys=True, indent=4, ensure_ascii=False)
    
    return scene,fathers

""" Print a model object in a geoJson file """
def printModelObject(theFile,tabs, i,name,verts,cells,father):
    tab = "    "
    print >> theFile, "-   ","id:", name
    print >> theFile, tab, "type:", "Feature"
    print >> theFile, tab, "geometry:" 
    print >> theFile, tab+tab, "type:", "Polygon"
    print >> theFile, tab+tab, "coordinates:" 
    print >> theFile, tab+tab+tab, "verts:"
    print >> theFile, tab+tab+tab+tab, AA(eval)(AA(vcode)(verts))
    print >> theFile, tab+tab+tab, "cells:"
    print >> theFile, tab+tab+tab+tab, cells
    print >> theFile, tab, "properties:"
    print >> theFile, tab+tab, "class:", "room"
    print >> theFile, tab+tab, "parent:", father
    print >> theFile, tab+tab, "son:", i
    print >> theFile, tab+tab, "description:", name
    print >> theFile, tab+tab, "tVector:", [0, 0, 0]
    print >> theFile, tab+tab, "rVector:", [0, 0, 0]

""" Print a mat object in a geoJson file """
def printMatObject(theFile,tabs, theMat):
    tab = "    "
    print >> theFile, tab+tab, "tVector:", theMat.T[-1,:-1].tolist()

""" Print a struct object in a geoJson file """
def printStructObject(theFile,tabs, i,name,box,father):
    tab = "    "
    print >> theFile, "-   ","id:", name
    print >> theFile, tab,"type:", "Feature"
    print >> theFile, tab,"geometry:" 
    print >> theFile, tab+tab, "type:", "Polygon"
    print >> theFile, tab+tab, "box:", box
    print >> theFile, tab,"properties:"
    print >> theFile, tab+tab, "class:", "assembly"
    print >> theFile, tab+tab, "parent:", father
    print >> theFile, tab+tab, "son:", i
    print >> theFile, tab+tab, "description:", name
    print >> theFile, tab+tab, "tVector:", [0, 0, 0]
    print >> theFile, tab+tab, "rVector:", [0, 0, 0]

""" Traverse a structure to print a geoJson file """
def printTraversal(theFile,CTM, stack, obj, scene=[], level=0, fathers=[],father=""):
   tabs = (4*level)*" "
   for i in range(len(obj)):
      if isinstance(obj[i],Model): 
         i,verts,cells = obj[i]
         if obj[i].__name__() == None:
            name = father+'-'+ str(id(obj[i]))
         else: 
            name = father+'-'+ str(obj[i].__name__())
         printModelObject(theFile,tabs, i,name,verts,cells,father)
         scene += [larApply(CTM)(obj[i])]
         fathers += [father]
      elif (isinstance(obj[i],tuple) or isinstance(obj[i],list)) and len(obj[i])==2:
         verts,cells = obj[i]
         name = father+'-'+ str(id(obj[i]))
         printModelObject(theFile,tabs, i,name,verts,cells,father)
         scene += [larApply(CTM)(obj[i])]
         fathers += [father]
      elif isinstance(obj[i],Mat): 
         printMatObject(theFile,tabs, obj[i])
         CTM = scipy.dot(CTM, obj[i])
      elif isinstance(obj[i], Struct):
         if obj[i].__name__() == None:
            name = father+'-'+ str(id(obj[i]))
         else: 
            name = father+'-'+ str(obj[i].__name__())
         box = obj[i].box
         printStructObject(theFile,tabs, i,name,box,father)
         stack.append(CTM) 
         level += 1
         fathers.append(name)
         printTraversal(theFile,CTM, stack, obj[i], scene, level, fathers,name)
         name = fathers.pop()
         level -= 1
         CTM = stack.pop()
   return scene,fathers


