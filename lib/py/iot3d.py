"""Module with automatic generation of simplified 3D buildings"""
import sys; sys.path.insert(0, 'lib/py/')
from architectural import *
import hijson

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
    V,boundaryEdges = hijson.structBoundaryModel(struct)   #+new
    boundaryPolygons = hijson.boundaryModel2polylines((V,boundaryEdges))  #+new
    
    p1V,p1FV,p1EV = struct2lar(struct)
    p1BE = AA(tuple)(AA(sorted)(boundaryEdges))
    p1IE = list(set(p1EV).difference(p1BE))
    
    #VIEW(STRUCT(MKPOLS((V,p1BE))))
    #VIEW(STRUCT(MKPOLS((V,p1IE))))

    printStructObject(theFile,tabs, 0,struct.name,struct.category,boundaryPolygons, "building")
    scene,fathers = printTraversal(theFile, CTM, stack, struct, [], 0, fathers,filename) 
    printWalls( theFile, p1V,p1BE,p1IE )
    theFile.close()
    """ exporting to JSON via YML a lar structure """
    file_handle = open(path+filename+".yml")
    my_dictionary = yaml.safe_load(file_handle)
    file_handle.close()
    with open(path+filename+".json", 'w') as outfile:
         json.dump(my_dictionary, outfile, sort_keys=True, indent=4, ensure_ascii=False)
    
    return scene,fathers

""" Print a model object in a geoJson file """
def printModelObject(theFile,tabs, i,name,category,verts,cells,father,tvect=[0,0,0]):
    tab = "    "
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
""" Print a model object in a geoJson file """
def printModelObject(theFile,tabs, i,name,category,verts,cells,father,tvect=[0,0,0]):
    tab = "    "
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
def printWalls( theFile, p1V,p1BE,p1IE ):
   
   for k,edge in enumerate(p1BE):
      tvectBE = VECTDIFF([[0,0,0], p1V[edge[0]]])
      edgeBE = [VECTDIFF([p1V[edge[0]],tvectBE]), VECTDIFF([p1V[edge[1]],tvectBE])]
      name = "e_wall_" + str(k)
      print "tvect=",tvectBE," edge=",edgeBE
      printWallObject(theFile, name,k, 'external_wall',edgeBE,'piano13D',tvectBE)
      
   for k,edge in enumerate(p1IE):
      tvectIE = VECTDIFF([[0,0,0], p1V[edge[0]]])
      edgeIE = [VECTDIFF([p1V[edge[0]],tvectIE]), VECTDIFF([p1V[edge[1]],tvectIE])]
      name = "i_wall_" + str(k)
      print "tvect=",tvectIE," edge=",edgeIE
      printWallObject(theFile, name,k, 'internal_wall',edgeIE,'piano13D',tvectIE)

""" Print a struct object in a geoJson file """
def printStructObject(theFile,tabs, i,name,category, boundaryPolyline,father,tvect=[0,0,0]):
    tab = "    "
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
    tab = "    "
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
def printTraversal(theFile,CTM, stack, obj, scene=[], level=0, fathers=[],father="",tvect=[0,0,0]):
   tabs = (4*level)*" "
   for i in range(len(obj)):
      if isinstance(obj[i],Model): 
         i,verts,cells = obj[i]
         if obj[i].__name__() == None:
            name = father+'-'+ str(id(obj[i]))
         else: 
            name = father+'-'+ str(obj[i].__name__())
         printModelObject(theFile,tabs, i,name,None,verts,cells,father)
         scene += [ rApply(CTM)(obj[i])]
         fathers += [father]
      elif (isinstance(obj[i],tuple) or isinstance(obj[i],list)) and len(obj[i])==2:
         verts,cells = obj[i]
         name = father+'-'+ str(id(obj[i]))
         printModelObject(theFile,tabs, i,name,None,verts,cells,father)
         scene += [larApply(CTM)(obj[i])]
         fathers += [father]
      elif isinstance(obj[i],Mat): 
         printMatObject(theFile,tabs, obj[i])
         CTM = scipy.dot(CTM,  obj[i])
      elif isinstance(obj[i], Struct):
         if obj[i].__name__() == None:
            name = father+'-'+ str(id(obj[i]))
         else: 
            name = father+'-'+ str(obj[i].__name__())
            
         box = obj[i].box
         tvect = box[0]
         V,boundaryEdges = hijson.structBoundaryModel(obj[i])   #+new
         #import pdb; pdb.set_trace()
         coordinates = hijson.boundaryModel2polylines((V,boundaryEdges))  #+new
         transfMat = array(t(*[-x for x in tvect]))
         
         def transformCycle(transfMat,cycle):
             affineCoordinates = [point+[1] for point in cycle]
             localCoords = (transfMat.dot(array(affineCoordinates).T)).T
             print localCoords
             localCoordinates = [[x,y,z] for x,y,z,w in localCoords.tolist()]
             return localCoordinates
             
         coordinates = [transformCycle(transfMat,cycle) for cycle in coordinates]
         category = obj[i].__category__()
         printStructObject(theFile,tabs, i,name,category,coordinates,father,tvect)
         stack.append(CTM) 
         level += 1
         fathers.append(name)
         printTraversal(theFile,CTM, stack, obj[i], scene, level, fathers,name, tvect)
         name = fathers.pop()
         level -= 1
         CTM = stack.pop()
   return scene,fathers


