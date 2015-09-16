""" sysml library """
from larlib import *
DRAW = COMP([VIEW,STRUCT,MKPOLS])


def boundaryOfChain(cells,facets):
   csrBoundaryMat = boundary(cells,facets)
   csrChain = zeros((len(cells),1))
   def boundaryOfChain0(chain):
      for cell in chain:  csrChain[cell,0]=1.0
      csrBoundaryChain = matrixProduct(csrBoundaryMat, csrChain)
      boundaryCells = [k for k,val in enumerate(csrBoundaryChain.tolist()) 
                     if val == [1.0]]
      return boundaryCells
   return boundaryOfChain0

""" Diagram initialization """
def assemblyDiagramInit (shape):
   def assemblyDiagram (quoteList):
      print "\n shape =",shape
      # shape and quoteList must be 3D, i.e. a python array with 3 indices
      assert (len(shape) == 3) and (len(quoteList) == 3)
      coordList = [list(cumsum([0]+pattern)) for pattern in quoteList]
      verts = CART(coordList)
      _,CV = larCuboids(shape)
      return verts,CV
   return assemblyDiagram

def lar2boundaryFaces(CV,FV):
   """ Boundary cells computation """
   return boundaryCells(CV,FV)

def lar2InteriorFaces(CV,FV):
   """ Boundary cells computation """
   boundarychain2D = boundaryCells(CV,FV)
   totalChain2D = range(len(FV))
   interiorCells = set(totalChain2D).difference(boundarychain2D)
   return interiorCells

""" Diagram scaling to given size """
def unitDiagram(diagram, size=[1,1,1]):
   V,CV = diagram
   print "\n shape =",shape
   # size must be a python array with 3 numbers
   assert (len(size) == 3) and (AND(AA(ISNUM)(size)) == True)
   V_ = array(V) / AA(float)(max(V))
   V = (V_ * size).tolist()
   diagram = V,CV
   return diagram



""" 3D window to viewport transformation """
def diagram2cellMatrix(diagram):
   def diagramToCellMatrix0(master,cell):
      wdw = min(diagram[0]) + max(diagram[0])         # window3D
      cV = [master[0][v] for v in master[1][cell]]
      vpt = min(cV) + max(cV)                      # viewport3D
      print "\n window3D =",wdw
      print "\n viewport3D =",vpt
      
      mat = zeros((4,4))
      mat[0,0] = (vpt[3]-vpt[0])/(wdw[3]-wdw[0])
      mat[0,3] = vpt[0] - mat[0,0]*wdw[0]
      mat[1,1] = (vpt[4]-vpt[1])/(wdw[4]-wdw[1])
      mat[1,3] = vpt[1] - mat[1,1]*wdw[1]
      mat[2,2] = (vpt[5]-vpt[2])/(wdw[5]-wdw[2])
      mat[2,3] = vpt[2] - mat[2,2]*wdw[2]
      mat[3,3] = 1
      print "\n mat =",mat
      return mat
   return diagramToCellMatrix0


def diagram2cell(diagram,master,cell):
   mat = diagram2cellMatrix(diagram)(master,cell)
   diagram =larApply(mat)(diagram)  
   (V1,CV1),(V2,CV2) = master,diagram
   n1,n2 = len(V1), len(V2)
   
   # identification of common vertices
   V, CV1, CV2, n12 = vertexSieve(master,diagram)
   commonRange = range(n1-n12, n1)
   newRange = range(n1,n1-n12+n2)
   
   # addition of incident vertices into the adjacents of theCell
   def checkInclusion(V,theCell,newRange):
      theVerts = [V[v] for v in theCell]
      theMin, theMax = min(theVerts), max(theVerts)
      theCell += [v for v in newRange if (
         theMin[0] <= V[v][0] and theMin[1] <= V[v][1] and theMin[2] <= V[v][2] 
         and 
         V[v][0] <= theMax[0] and V[v][1] <= theMax[1] and V[v][2] <= theMax[2] 
         )]
      return theCell
   
   # addition of new vertices into the adjacents of cell c
   CV1 = [checkInclusion(V,c,newRange) 
         if set(c).intersection(commonRange) != set() else c
          for c in CV1]
   
   # masterBoundaryFaces = boundaryOfChain(CV,FV)([cell])
   # diagramBoundaryFaces = lar2boundaryFaces(CV,FV)
   CV = [c for k,c in enumerate(CV1) if k != cell] + CV2
   
   master = V, CV
   return master

""" Exterior space of a block diagram """
def exteriorCells(diagram):
   V,CV = diagram
   minVert, maxVert = min(V), max(V)
   d = len(V[0])
   outchain = [[] for k in range(2*d)]
   for k,v in enumerate(V):
      for h in range(d):
         if v[h] == minVert[h]: outchain[h] += [k]
         if v[h] == maxVert[h]: outchain[h+d] += [k]
   return outchain

""" Place vertices of two LAR models in a common space """
from collections import defaultdict, OrderedDict

def vertexSieve(model1, model2):
   V1,CV1 = model1; V2,CV2 = model2
   n = len(V1); m = len(V2)
   def shift(CV, n): 
      return [[v+n for v in cell]for cell in CV]
   CV2 = shift(CV2,n)


   vdict1 = defaultdict(list)
   for k,v in enumerate(V1): vdict1[vcode(v)].append(k) 
   vdict2 = defaultdict(list)
   for k,v in enumerate(V2): vdict2[vcode(v)].append(k+n) 
   
   vertdict = defaultdict(list)
   for point in vdict1.keys(): vertdict[point] += vdict1[point]
   for point in vdict2.keys(): vertdict[point] += vdict2[point]


   case1, case12, case2 = [],[],[]
   for item in vertdict.items():
      key,val = item
      if len(val)==2:  case12 += [item]
      elif val[0] < n: case1 += [item]
      else: case2 += [item]
   n1 = len(case1); n2 = len(case12); n3 = len(case2)


   invertedindex = list(0 for k in range(n+m))
   for k,item in enumerate(case1):
      invertedindex[item[1][0]] = k
   for k,item in enumerate(case12):
      invertedindex[item[1][0]] = k+n1
      invertedindex[item[1][1]] = k+n1
   for k,item in enumerate(case2):
      invertedindex[item[1][0]] = k+n1+n2


   V = [eval(p[0]) for p in case1] + [eval(p[0]) for p in case12] + [eval(
            p[0]) for p in case2]
   CV1 = [sorted([invertedindex[v] for v in cell]) for cell in CV1]
   CV2 = [sorted([invertedindex[v] for v in cell]) for cell in CV2]
   return V, CV1, CV2, len(case12)


