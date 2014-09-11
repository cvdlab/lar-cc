""" Module for Boolean ops with LAR """
from pyplasm import *
from scipy import *
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from lar2psm import *
from simplexn import *
from larcc import *
from largrid import *
from myfont import *
from mapper import *

DEBUG = True
from matrix import *
from splitcell import *
""" TODO: use package Decimal (http://docs.python.org/2/library/decimal.html) """
global PRECISION
PRECISION = 4.95

def verySmall(number): return abs(number) < 10**-(PRECISION/1.15)

def prepKey (args): return "["+", ".join(args)+"]"

def fixedPrec(value):
   out = round(value*10**(PRECISION*1.1))/10**(PRECISION*1.1)
   if out == -0.0: out = 0.0
   return str(out)
   
def vcode (vect): 
   """
   To generate a string representation of a number array.
   Used to generate the vertex keys in PointSet dictionary, and other similar operations.
   """
   return prepKey(AA(fixedPrec)(vect))

""" First step of Boolean Algorithm """
from collections import defaultdict, OrderedDict

""" TODO: change defaultdict to OrderedDefaultdict """

class OrderedDefaultdict(collections.OrderedDict):
   def __init__(self, *args, **kwargs):
      if not args:
         self.default_factory = None
      else:
         if not (args[0] is None or callable(args[0])):
            raise TypeError('first argument must be callable or None')
         self.default_factory = args[0]
         args = args[1:]
      super(OrderedDefaultdict, self).__init__(*args, **kwargs)

   def __missing__ (self, key):
      if self.default_factory is None:
         raise KeyError(key)
      self[key] = default = self.default_factory()
      return default

   def __reduce__(self):  # optional, for pickle support
      args = (self.default_factory,) if self.default_factory else tuple()
      return self.__class__, args, None, None, self.iteritems()


def vertexSieve(model1, model2):
   from lar2psm import larModelBreak
   V1,CV1 = larModelBreak(model1) 
   V2,CV2 = larModelBreak(model2)
   n = len(V1); m = len(V2)
   def shift(CV, n): 
      return [[v+n for v in cell] for cell in CV]
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


def covering(model1,model2,dim=2,emptyCellNumber=1):
   V, CV1, CV2, n12 = vertexSieve(model1,model2)
   _,EEV1 = larFacets((V,CV1),dim,emptyCellNumber)
   _,EEV2 = larFacets((V,CV2),dim,emptyCellNumber)
   if emptyCellNumber !=0: CV1 = CV1[:-emptyCellNumber]
   if emptyCellNumber !=0: CV2 = CV2[:-emptyCellNumber]
   VV = AA(LIST)(range(len(V)))
   return V,[VV,EEV1,EEV2,CV1,CV2],n12


""" Characteristic matrix transposition """
def invertRelation(dim,CV):
   inverse = [[] for k in range(dim)]
   for k,cell in enumerate(CV):
      for v in cell:
         inverse[v] += [k]
   return inverse





""" Cell splitting in two cells """
def cellSplitting(face,cell,covector,V,EEV,CV):
   plane = COVECTOR([V[v] for v in EEV[face]])
   theCell = [V[v] for v in CV[cell]]
   [below,equal,above] = SPLITCELL(plane,theCell)
   cell1 = AA(vcode)(below)
   cell2 = AA(vcode)(above)
   return cell1,cell2

""" Init face-cell and cell-face dictionaries """
def initTasks(tasks):
   dict_fc = defaultdict(list)
   dict_cf = defaultdict(list)
   for task in tasks:
      face,cell,covector = task
      dict_fc[face] += [(cell,covector)] 
      dict_cf[cell] += [(face,covector)] 
   return dict_fc,dict_cf


""" Updating the vertex set of split cells """

""" Computation of bits of split cells """
def testingSubspace(V,covector):
   def testingSubspace0(vcell):
      inout = SIGN(sum([INNERPROD([V[v]+[1.],covector]) for v in vcell]))
      return inout
   return testingSubspace0
   
def cuttingTest(covector,polytope,V):
   signs = [INNERPROD([covector, V[v]+[1.]]) for v in polytope]
   signs = eval(vcode(signs))
   return any([value<-0.001 for value in signs]) and any([value>0.001 for value in signs])


def tangentTest(face,polytope,V,BC):
   faceVerts = BC[face]
   cellVerts = polytope
   commonVerts = list(set(faceVerts).intersection(cellVerts))
   if commonVerts != []:
      v0 = commonVerts[0] # v0 = common vertex (TODO more general)
      transformMat = mat([DIFF([V[v],V[v0]]) for v in cellVerts if v != v0]).T.I
      vects = (transformMat * (mat([DIFF([V[v],V[v0]]) for v in faceVerts 
               if v != v0]).T)).T.tolist()
      if all([all([x>=-0.0001 for x in list(vect)]) for vect in vects]): 
         return True
   else: return False

from collections import defaultdict

def splitCellsCreateVertices(vertdict,dict_fc,dict_cf,V,BC,CV,VC,lenBC1):
   DEBUG = False
   splitBoundaryFacets = []; splittingCovectors = defaultdict(list)
   CVbits = [[-1,-1] for k in range(len(CV))] 
   nverts = len(V); cellPairs = []; twoCellIndices = []; 
   while any([tasks != [] for face,tasks in dict_fc.items()]) : 
      for face,tasks in dict_fc.items():
         for task in tasks:
            cell,covector = task
            vcell = CV[cell]
            if (cell,vcell,face) == (29, [24, 1, 4], 3): break
            print "\n1> cell,vcell,face,covector =",cell,vcell,face,covector

            cell1,cell2 = cellSplitting(face,cell,covector,V,BC,CV)
            if cuttingTest(covector,vcell,V):
               if cell1 == [] or cell2 == []:
                  pass
               else:
                  adjCells = adjacencyQuery(V,CV)(cell)
                                    
                  vcell1 = []
                  for k in cell1:
                     if vertdict[k]==[]: 
                        vertdict[k] += [nverts]
                        V += [eval(k)]
                        nverts += 1
                     vcell1 += [vertdict[k]]
                  
                  vcell1 = CAT(vcell1)
                  vcell2 = CAT([vertdict[k] for k in cell2])                     
                                             
                  V,CV,CVbits, dict_cf, dict_fc,twoCells = splittingControl(
                     face,cell,covector,vcell,vcell1,vcell2, dict_fc,dict_cf,V,BC,CV,VC,
                     CVbits,lenBC1,splitBoundaryFacets,splittingCovectors)
                  if twoCells[0] != twoCells[1]:

                     print "2> cell,adjCells =",cell,adjCells
                     
                     for adjCell in adjCells:
                        if cuttingTest(covector,CV[adjCell],V):
                           dict_fc[face] += [(adjCell,covector)]                             
                           dict_cf[adjCell] += [(face,covector)] 
                           print "2.1> face,dict_fc[face] =",face,dict_fc[face]
                           print "2.2> adjCell,dict_cf[adjCell] =",adjCell,dict_cf[adjCell]

                        cellPairs += [[vcell1, vcell2]]
                        twoCellIndices += [[twoCells]]
                                    
               if DEBUG: showSplitting("c",twoCells[1],V,cellPairs,BC,CV)

            elif tangentTest(face,vcell,V,BC):                             
               newFacet = [ v for v in vcell if 
                  verySmall(INNERPROD([covector,V[v]+[1.0]])) ]
               splitBoundaryFacets += [newFacet]
               splittingCovectors[cell] += [(face,covector,newFacet)]
               
               dict_fc[face].remove((cell,covector))   # remove the split cell
               dict_cf[cell].remove((face,covector))   # remove the splitting face

            else: 
               dict_fc[face].remove((cell,covector))   # remove the split cell
               # dict_cf[cell].remove((face,covector))   # remove the splitting face
            if DEBUG: showSplitting("b",cell,V,cellPairs,BC,CV)
         if DEBUG: showSplitting("a",cell,V,cellPairs,BC,CV)
   splitBoundaryFacets = sorted(list(AA(list)(set(AA(tuple)(AA(sorted)(splitBoundaryFacets))))))
   return CVbits,cellPairs,twoCellIndices,splitBoundaryFacets,splittingCovectors

""" Managing the splitting dictionaries """
def splittingControl(face,cell,covector,vcell,vcell1,vcell2,
      dict_fc,dict_cf,V,BC,CV,VC,CVbits,lenBC1,splitBoundaryFacets,splittingCovectors):

   boundaryFacet = BC[face]
   translVector = V[boundaryFacet[0]]
   tcovector = [cv+tv*covector[-1] for (cv,tv) in zip(
               covector[:-1],translVector) ]+[0.0]

   c1,c2 = cell,cell
   if not haltingSplitTest(face,cell,vcell,vcell1,vcell2,boundaryFacet,
                        translVector,tcovector,covector,
                        V,splitBoundaryFacets,splittingCovectors) :
                        
      print "1.1> cell,vcell,face,covector =",cell,vcell,face,covector

      # only one facet covector crossing the cell
      cellVerts = CV[cell]
      CV[cell] = vcell1
      CV += [vcell2]
      CVbits += [list(copy(CVbits[cell]))]
      c1,c2 = cell,len(CV)-1
      
      newFacet = list(set(vcell1).intersection(vcell2))
      splitBoundaryFacets += [newFacet]  ## CAUTION: to verify
      splittingCovectors[c1] += [(face,covector,newFacet)]
      splittingCovectors[c2] = splittingCovectors[c1]
      
      print "1.1.1> c1,c2,CVbits[c1],CVbits[c2] =",c1,c2,CVbits[c1],CVbits[c2]
      
      dict_fc[face].remove((cell,covector))  # remove the split cell
      dict_cf[cell].remove((face,covector))  # remove the splitting face
            
      # more than one facet covectors crossing the cell
      alist1,alist2 = list(),list()
      for aface,covector in dict_cf[cell]:
         if cuttingTest(covector,CV[cell],V):
      
            # for each facet crossing the cell
            # compute the intersection between the facet and the cell
            faceVerts = BC[aface]
            commonVerts = list(set(faceVerts).intersection(cellVerts))
            
            # and attribute the intersection to the split subcells
            if set(vcell1).intersection(commonVerts) != set():
               alist1.append((aface,covector))
            else: dict_fc[aface].remove((cell,covector)) 
                  
            if set(vcell2).intersection(commonVerts) != set():
               alist2.append((aface,covector))
               dict_fc[aface] += [(len(CV)-1,covector)]
         
            print "1.1.1.1> aface,dict_fc[aface] =",aface,dict_fc[aface]
         
      dict_cf[cell] = alist1  
      dict_cf[len(CV)-1] = alist2
      
      print "1.1.2> cell,dict_cf[cell] =",cell,dict_cf[cell]
      print "1.1.3> len(CV)-1,dict_cf[len(CV)-1] =",len(CV)-1,dict_cf[len(CV)-1]
      
   else:
   
      print "1.2> cell,vcell,face,covector =",cell,vcell,face,covector
      
      dict_fc[face].remove((cell,covector))  # remove the split cell
      dict_cf[cell].remove((face,covector))  # remove the splitting face   
      
   return V,CV,CVbits, dict_cf, dict_fc,[c1,c2]

""" Test for split halting along a boundary facet """
def haltingSplitTest(face,cell,vcell,vcell1,vcell2,boundaryFacet,translVector,tcovector,covector,
                  V,splitBoundaryFacets,splittingCovectors):
   newFacet = list(set(vcell1).intersection(vcell2))
   
   # translation 
   newFacet = [ eval(vcode(VECTDIFF([V[v],translVector]))) for v in newFacet ]
   boundaryFacet = [ eval(vcode(VECTDIFF([V[v],translVector]))) for v in boundaryFacet ]
   
   # linear transformation: newFacet -> standard (d-1)-simplex
   transformMat = mat( boundaryFacet[1:] + [tcovector[:-1]] ).T.I
   
   # transformation in the subspace x_d = 0
   newFacet = AA(COMP([eval,vcode]))((transformMat * (mat(newFacet).T)).T.tolist())
   boundaryFacet = AA(COMP([eval,vcode]))((transformMat * (mat(boundaryFacet).T)).T.tolist())
   
   # projection in E^{d-1} space and Boolean test
   newFacet = MKPOL([ AA(lambda v: v[:-1])(newFacet), [range(1,len(newFacet)+1)], None ])
   boundaryFacet = MKPOL([ AA(lambda v: v[:-1])(boundaryFacet), [range(1,len(boundaryFacet)+1)], None ])
   verts,cells,pols = UKPOL(INTERSECTION([newFacet,boundaryFacet]))
   
   if verts == []: return True
   else: return False

""" Computing the adjacent cells of a given cell """
def adjacencyQuery (V,CV):
   dim = len(V[0])
   def adjacencyQuery0 (cell):
      nverts = len(CV[cell])
      csrCV =  csrCreate(CV)
      csrAdj = matrixProduct(csrCV,csrTranspose(csrCV))
      cellAdjacencies = csrAdj.indices[csrAdj.indptr[cell]:csrAdj.indptr[cell+1]]
      return [acell for acell in cellAdjacencies if dim <= csrAdj[cell,acell] < nverts]
   return adjacencyQuery0

""" Show the process of CDC splitting """
def showSplitting(step,theCell,V,cellPairs,BC,CV):
   VV = AA(LIST)(range(len(V)))
   boundaries = COLOR(RED)(SKEL_1(STRUCT(MKPOLS((V,BC)))))
   submodel = COLOR(CYAN)(STRUCT([ SKEL_1(STRUCT(MKPOLS((V,CV)))), boundaries ]))
   if cellPairs != []:
      cells1,cells2 = TRANS(cellPairs)
      out = [COLOR(WHITE)(MKPOL([V,[[v+1 for v in cell] for cell in cells1],None])), 
            COLOR(MAGENTA)(MKPOL([V,[[v+1 for v in cell] for cell in cells2],None]))]
      VIEW(STRUCT([ STRUCT(out),larModelNumbering(V,[VV,BC,CV],submodel,2), 
         S([1,2])([0.1,0.1])(TEXT(str(theCell)+step)) ]))
   else:
      VIEW(STRUCT([ larModelNumbering(V,[VV,BC,CV],submodel,2),
         S([1,2])([0.1,0.1])(TEXT(str(theCell)+step)) ]))

""" Boundary triangulation of a convex hull """
"""
def qhullBoundary(V):
   dim = len(V[0])
   triangulation = Delaunay(array(V))
   CV = triangulation.simplices.tolist()
   Ad = triangulation.neighbors.tolist()
   wingedRep = zip(CV,Ad)
   boundaryCofaces = [simplex for simplex in wingedRep if any([ad==-1 for ad in simplex[1]])]
   wingedPairs = [zip(*coface) for coface in boundaryCofaces]
   out = [[v for v,ad in pairs if ad!=-1] for pairs in wingedPairs]
   return sorted(out)
"""
from scipy.spatial import ConvexHull
def qhullBoundary(V):
   points = array(V)
   hull = ConvexHull(points)
   out = hull.simplices.tolist()
   return sorted(out)
   
if __name__=="__main__":
   BV = qhullBoundary(V)
   VIEW(STRUCT(MKPOLS((V,BV))))

""" Extracting a $(d-1)$-basis of SCDC """
def larConvexFacets (V,CV):
   dim = len(V[0])
   model = V,CV
   V,FV = larFacets(model,dim)
   FV = sorted(FV + qhullBoundary(V))
   return FV
   
if __name__=="__main__":
   V = [[0.0,10.0],[0.0,0.0],[10.0,10.0],[10.0,0.0],[12.5,2.5],[2.5,2.5],[2.5,12.5],
       [12.5,12.5],[10.0,2.5],[2.5,10.0]]
   CV = [[0,1,5],[9,0,5],[9,0,6],[1,3,5],[8,4,3],[8,5,3],[2,4,7],[2,6,7],
        [8,2,5],[8,4,2],[9,2,6],[9,2,5]]
   VV = AA(LIST)(range(len(V)))
   FV = larConvexFacets (V,CV)
   submodel = SKEL_1(STRUCT(MKPOLS((V,CV))))
   VIEW(larModelNumbering(V,[VV,FV,CV],submodel,4))

""" Traversing a Boolean argument within the CDC """
def booleanChainTraverse(h,cell,V,CV,CVbits,value):
   adjCells = adjacencyQuery(V,CV)(cell)
   for adjCell in adjCells: 
      if CVbits[adjCell][h] == -1:
         CVbits[adjCell][h] = value
         CVbits = booleanChainTraverse(h,adjCell,V,CV,CVbits,value)
   return CVbits

""" Boolean fragmentation and classification of CDC """

def booleanChains(arg1,arg2):
   (V1,basis1), (V2,basis2) = arg1,arg2
   model1, model2 = (V1,basis1[-1]), (V2,basis2[-1])
   V,[VV,_,_,CV1,CV2],n12 = covering(model1,model2,2,0)
   CV = sorted(AA(sorted)(Delaunay(array(V)).simplices))
   vertdict = defaultdict(list)
   for k,v in enumerate(V): vertdict[vcode(v)] += [k]
   
   BC1 = signedCellularBoundaryCells(V1,basis1)
   BC2 = signedCellularBoundaryCells(V2,basis2)
   n_bf1,n_bf2 = len(BC1),len(BC2)
   BC = [[ vertdict[vcode(V1[v])][0] for v in cell] for cell in BC1] + [ 
         [ vertdict[vcode(V2[v])][0] for v in cell] for cell in BC2]
   BV = list(set(CAT([v for v in BC])))
   VV = AA(LIST)(range(len(V)))
   
   if DEBUG: 
      """ Input and CDC visualisation """
      dim = len(V[0])
      if dim == 2:
          submodel1 = mkSignedEdges((V1,BC1))
          submodel2 = mkSignedEdges((V2,BC2))
          VIEW(STRUCT([submodel1,submodel2]))
      submodel = SKEL_1(STRUCT(MKPOLS((V,CV))))
      VIEW(larModelNumbering(V,[VV,[],CV],submodel,4))
      submodel = STRUCT([SKEL_1(STRUCT(MKPOLS((V,CV)))), COLOR(RED)(STRUCT(MKPOLS((V,BC))))])
      VIEW(larModelNumbering(V,[VV,BC,CV],submodel,4))
      
      
   """ New implementation of splitting dictionaries """
   VC = invertRelation(len(V),CV)
   
   covectors = []
   for faceVerts in BC:
      points = [V[v] for v in faceVerts]
      """
      dim = len(points[0])
      theMat = Matrix( [(dim+1)*[1.]] + [p+[1.] for p in points] )
      covector1 = [(-1)**(col)*theMat.minor(0,col).determinant() 
                     for col in range(dim+1)]
      """
      covector = COVECTOR(points)
      covector2 = covector[1:]+[covector[0]] 
      print "covector =",covector2
      covectors += [covector2]
   
   """ to compute a single d-cell associated to (face,covector) """
   def covectorCell(face,faceVerts,covector,CV,VC):
      incidentCells = VC[faceVerts[0]]
      for cell in incidentCells:
         cellVerts = CV[cell]
         v0 = list(set(faceVerts).intersection(cellVerts))[0] # v0 = common vertex
         transformMat = mat([DIFF([V[v],V[v0]]) for v in cellVerts if v != v0]).T.I
         vects = (transformMat * (mat([DIFF([V[v],V[v0]]) for v in faceVerts 
                  if v != v0]).T)).T.tolist()
         if any([all([x>=-0.0001 for x in list(vect)]) for vect in vects]): 
            return [face,cell,covector]
      print "error: found no face,cell,covector","\n"
   
   """ Initialization of splitting dictionaries """
   tasks = []
   for face,covector in zip(range(len(BC)),covectors):
      tasks += [covectorCell(face,BC[face],covector,CV,VC)]
   
   dict_fc,dict_cf = initTasks(tasks)
   print "\ndict_fc =",dict_fc
   print "dict_cf =",dict_cf,"\n"
   
   
   
   CVbits,cellPairs,twoCellIndices,splitBoundaryFacets,splittingCovectors = \
      splitCellsCreateVertices( vertdict,dict_fc,dict_cf,V,BC,CV,VC,len(BC1) )
   showSplitting("z",len(CV),V,cellPairs,BC,CV)
   
   """ Numerical instability of vertices curation """
   dim = len(V[0])
   if dim == 2:
       x,y = TRANS(V)
       tree = scipy.spatial.KDTree(zip(array(x).ravel(), array(y).ravel()))
   elif dim == 3:
       x,y,z = TRANS(V)
       tree = scipy.spatial.KDTree(zip(array(x).ravel(), array(y).ravel(), array(z).ravel()))
   closestVertexPairs = AA(list)(tree.query(tree.data,2)[1])
   distances = sorted([[VECTNORM(VECTDIFF([V[v],V[w]])),v,w] for v,w in closestVertexPairs])
   coincidentVertexPairs = [[v,w] for k,(dist,v,w) in enumerate(distances) if dist < 10**-PRECISION]
   
   # remove w from CV (v <- w)
   if coincidentVertexPairs != []:
      coincidentVertexPairs = list(set(AA(tuple)(AA(sorted)(coincidentVertexPairs))))
      toChange = TRANS(coincidentVertexPairs)[1]
      mapping = dict(AA(REVERSE)(coincidentVertexPairs))
      CV_ = [[v  if v not in toChange else mapping[v] for v in cell] for cell in CV]
      VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,larConvexFacets (V,CV_)))))
      CV = CV_
   

   print "\ndict_fc =",dict_fc
   print "dict_cf =",dict_cf,"\n"

   """ Building a dictionary of SCDC $(d-1)$-cells """
   def facetBasisDict(model):
      V,CV = model
      FV = larConvexFacets (V,CV)
      values = range(len(FV))
      keys = AA(tuple)(FV)
      dict_facets = dict(zip(keys,values))
      return dict_facets
      
   """ Searching for the split boundary facets in the dictionary """
   if __name__=="__main__":
      model = V,CV
      dict_facets = facetBasisDict(model)
      for cell in splitBoundaryFacets: 
         if cell in dict_facets:
            print dict_facets[cell]
         else: print cell
   
   
   dict_facets = facetBasisDict((V,CV))
   for cell in AA(tuple)(splitBoundaryFacets): 
      if cell in dict_facets:
         print dict_facets[cell]
      else: print cell
      
   VV = AA(LIST)(range(len(V)))  
   submodel = STRUCT(MKPOLS((V,larConvexFacets (V,CV))))
   VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,larConvexFacets (V,CV)))))
   
   return V,CV,BC,CVbits,vertdict,dict_facets,splittingCovectors,n_bf1,n_bf2

""" Extraction of LAR reps of common Boolean status """
def larBooleanPartition(CVbits,CV):
   ordCV = sorted(zip(CVbits,CV))
   out = defaultdict(list)
   for status,cell in ordCV:
      out[tuple(status)] += [cell]
   return out

