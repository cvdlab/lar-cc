"""Module with automatic generation of simplified 3D buildings"""
from larlib import *


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
   #import pdb; pdb.set_trace()

   tabs = (3*1)*" "

   import hijson
   V,boundaryEdges = hijson.structBoundaryModel(struct)   #+new
   boundaryPolygons = hijson.boundaryModel2polylines((V,boundaryEdges))  #+new
   
   printStructObject(theFile,tabs, 0,struct.name,struct.category,boundaryPolygons, "building")
   printTraversal(theFile, struct, 0, fathers,filename,[0,0,0]) 

   theFile.close()
   """ exporting to JSON via YML a lar structure """
   file_handle = open(path+filename+".yml")
   my_dictionary = yaml.safe_load(file_handle)
   file_handle.close()
   with open(path+filename+".json", 'w') as outfile:
       json.dump(my_dictionary, outfile, sort_keys=True, indent=4, ensure_ascii=False)
   #filejson = open("file1.json", 'r')
   #dicta = json.loads(filejson.read())
   

""" Print a model object in a geoJson file """
def printModelObject(theFile,tabs, i,name,category,verts,cells,father,tvect=[0,0,0]):
   tab = "  "
   print >> theFile, "-   ","id:", name
   print >> theFile, tab, "type:", "Feature"
   print >> theFile, tab, "geometry:" 
   print >> theFile, tab+tab, "type:", "Polygon"
   print >> theFile, tab+tab, "coordinates:" 
   print >> theFile, tab+tab+tab, AA(eval)(AA(vcode)(verts))
   print >> theFile, tab, "properties:"
   print >> theFile, tab+tab, "class:", category
   print >> theFile, tab+tab, "parent:", father
   print >> theFile, tab+tab, "son:", i
   print >> theFile, tab+tab, "description:", name
   print >> theFile, tab+tab, "tVector:", [tvect[0],tvect[1],tvect[2]]
   print >> theFile, tab+tab, "rVector:", [0, 0, 0]

""" Print a mat object in a geoJson file """
def printMatObject(theFile,tabs, theMat):
   tab = "  "
   print >> theFile, tab+tab, "tVector:", theMat.T[-1,:-1].tolist()

""" Print a struct object in a geoJson file """
def printStructObject(theFile,tabs, i,name,category, boundaryPolyline,father,tvect=[0,0,0]):
   tab = "  "
   print >> theFile, "-   ","id:", name
   print >> theFile, tab,"type:", "Feature"
   print >> theFile, tab,"geometry:" 
   print >> theFile, tab+tab, "type:", "Polygon"
   print >> theFile, tab+tab, "coordinates:", boundaryPolyline
   print >> theFile, tab,"properties:"
   print "category=", category
   print >> theFile, tab+tab, "class:", category
   print >> theFile, tab+tab, "parent:", father
   print >> theFile, tab+tab, "son:", i
   print >> theFile, tab+tab, "description:", name
   print >> theFile, tab+tab, "tVector:", [tvect[0],tvect[1],tvect[2]]
   print >> theFile, tab+tab, "rVector:", [0, 0, 0]
   
""" Print a struct object in a geoJson file """
def printWallObject(theFile, name,i,category, edge,father,tvect=[0,0,0]):
   tab = "  "
   print >> theFile, "-   ","id:", name
   print >> theFile, tab,"type:", "Feature"
   print >> theFile, tab,"geometry:" 
   print >> theFile, tab+tab, "type:", "LineString"
   print >> theFile, tab+tab, "coordinates:", edge
   print >> theFile, tab,"properties:"
   print "category=", category
   print >> theFile, tab+tab, "class:", category
   print >> theFile, tab+tab, "parent:", father
   print >> theFile, tab+tab, "son:", i
   print >> theFile, tab+tab, "description:", name
   print >> theFile, tab+tab, "tVector:", [tvect[0],tvect[1],tvect[2]]
   print >> theFile, tab+tab, "rVector:", [0, 0, 0]


""" Traverse a structure to print a geoJson file """
def old_printTraversal(theFile,obj, level=0, fathers=[],father="",tvect=[0,0,0]):
   import hijson
   tabs = (4*level)*" "
   for i in range(len(obj)):
     if isinstance(obj[i],Model): 
       i,verts,cells = obj[i]
       if obj[i].__name__() == None:
         name = father+'-'+ str(id(obj[i]))
       else: 
         name = father+'-'+ str(obj[i].__name__())
       printModelObject(theFile,tabs, i,name,None,verts,cells,father)
       fathers += [father]
     elif (isinstance(obj[i],tuple) or isinstance(obj[i],list)) and len(obj[i])==2:
       verts,cells = obj[i]
       name = father+'-'+ str(id(obj[i]))
       printModelObject(theFile,tabs, i,name,None,verts,cells,father)
       fathers += [father]
     elif isinstance(obj[i],Mat): 
       printMatObject(theFile,tabs, obj[i])
     elif isinstance(obj[i], Struct):
       if obj[i].__name__() == None:
         name = father+'-'+ str(id(obj[i]))
       else: 
         name = father+'-'+ str(obj[i].__name__())
         
       box = obj[i].box
       tvect = box[0]
       
       if obj[i].category in ["internal_wall", "external_wall", "corridor_wall"]:
          edge = evalStruct(obj[i])[0][0]
          category = obj[i].__category__()
          printWallObject(theFile, name,i,category, edge,father,tvect=[0,0,0])
       else:      
          V,boundaryEdges = hijson.structBoundaryModel(obj[i])   
          coordinates = hijson.boundaryModel2polylines((V,boundaryEdges))
          transfMat = array(t(*[-x for x in tvect]))
          
          def transformCycle(transfMat,cycle):
             affineCoordinates = [point+[1] for point in cycle]
             localCoords = (transfMat.dot(array(affineCoordinates).T)).T
             localCoordinates = [[x,y,z] for x,y,z,w in localCoords.tolist()]
             return localCoordinates
             
          coordinates = [transformCycle(transfMat,cycle) for cycle in coordinates]
          category = obj[i].__category__()
          printStructObject(theFile,tabs, i,name,category,coordinates,father,tvect)
          
       level += 1
       fathers.append(name)
       printTraversal(theFile,obj[i], level, fathers,name, tvect)
       name = fathers.pop()
       level -= 1
   return fathers



###

wall_categories = ["internal_wall", "external_wall", "corridor_wall"]

room_categories = ["room"]

indent = 2 * " "

def get_id(obj, father_id):
   return father_id + '_' + obj.name

def is_struct(obj):
   return isinstance(obj, Struct)
   
def is_wall(obj):
   return obj.category in wall_categories

def is_room(obj):
   return obj.category in room_categories

def is_model(obj):
   return isinstance(obj, Model)

def is_mat(obj):
   return isinstance(obj,Mat)

def is_list(obj):
   return isinstance(obj, tuple) or isinstance(obj, list)

def print_wall(out, obj, parent_id):
   obj_id = get_id(obj, parent_id)
   edge = evalStruct(obj)[0][0]
   category = obj.category
   tvect = obj.box[0]
   print >> out, "-   ", "id:", obj_id
   print >> out, indent, "type:", "Feature"
   print >> out, indent, "geometry:" 
   print >> out, indent+indent, "type:", "LineString"
   print >> out, indent+indent, "coordinates:", edge
   print >> out, indent, "properties:"
   print "category=", category
   print >> out, indent+indent, "class:", category
   print >> out, indent+indent, "parent:", parent_id
   print >> out, indent+indent, "description:", ""
   print >> out, indent+indent, "tVector:", [tvect[0],tvect[1],tvect[2]]
   print >> out, indent+indent, "rVector:", [0, 0, 0]  


def transformCycle(transfMat,cycle):
   affineCoordinates = [point+[1] for point in cycle]
   localCoords = (transfMat.dot(array(affineCoordinates).T)).T
   localCoordinates = [[x,y,z] for x,y,z,w in localCoords.tolist()]
   return localCoordinates

def print_struct(out, obj, parent_id):
   tvect = obj.box[0]
   import hijson
   V,boundaryEdges = hijson.structBoundaryModel(obj)   
   coordinates = hijson.boundaryModel2polylines((V,boundaryEdges))
   transfMat = array(t(*[-x for x in tvect]))   
   coordinates = [transformCycle(transfMat,cycle) for cycle in coordinates]
   category = obj.category
   
   print >> out, "-   ","id:", name
   print >> out, indent,"type:", "Feature"
   print >> out, indent,"geometry:" 
   print >> out, indent+indent, "type:", "Polygon"
   print >> out, indent+indent, "coordinates:", coordinates
   print >> out, indent,"properties:"
   print "category=", category
   print >> out, indent+indent, "class:", category
   print >> out, indent+indent, "parent:", parent_id
   print >> out, indent+indent, "description:", ""
   print >> out, indent+indent, "tVector:", [tvect[0],tvect[1],tvect[2]]
   print >> out, indent+indent, "rVector:", [0, 0, 0]

def print_model(out, obj, parent_id):
   obj_id = get_id(obj, parent_id)
   category = obj.category
   verts,cells = obj
   tvect = obj.box[0]
   print >> out, "-   ", "id:", obj_id
   print >> out, indent, "type:", "Feature"
   print >> out, indent, "geometry:" 
   print >> out, indent+indent, "type:", "Polygon"
   print >> out, indent+indent, "coordinates:", AA(eval)(AA(vcode)(verts))
   print >> out, indent, "properties:"
   print >> out, indent+indent, "class:", category
   print >> out, indent+indent, "parent:", parent_id
   print >> out, indent+indent, "description:", ""
   print >> out, indent+indent, "tVector:", [tvect[0],tvect[1],tvect[2]]
   print >> out, indent+indent, "rVector:", [0, 0, 0]

def print_mat(out, obj, parent_id):
   print >> out, indent+indent, "tVector:", obj.T[-1,:-1].tolist()


def print_traversal(out, obj, parent_id='', level=0): 

   if is_struct(obj):

      if is_wall(obj):
         print_wall(out, obj, parent_id)
         return

      print_struct(out, obj, parent_id)

      parent_id = get_id(obj, parent_id)
      print_traversal(out, obj.body, parent_id)
      return

   if is_model(obj):
      print_model(out, obj, parent_id)
      return

   if is_list(obj):
      parent_id = get_id(obj, parent_id)      
      print_traversal(out, obj[0], parent_id)
      
      if len(obj) > 1:
         print_traversal(out, obj[1:], parent_id)
   
      return
   
   if is_mat(obj):
      print_mat(out, obj, parent_id)
      return

def printStruct2GeoJson(path, struct):
   if struct.__name__() == None:
      filename = str(id(struct))
   else: 
      filename = struct.__name__()
   out = open(path+filename+".yml", "w")
   print >> out, "---"
   print "filename =", path+filename+".yml"
   dim = checkStruct(struct.body)
   CTM, stack = scipy.identity(dim+1), []
   fathers = [filename]
   #import pdb; pdb.set_trace()

   import hijson
   V,boundaryEdges = hijson.structBoundaryModel(struct)   #+new
   boundaryPolygons = hijson.boundaryModel2polylines((V,boundaryEdges))  #+new
   
   print_traversal(out, struct) 

   out.close()
   


