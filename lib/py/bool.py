""" Module for Boolean ops with LAR """
from matrix import *
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

""" TODO: use package Decimal (http://docs.python.org/2/library/decimal.html) """
PRECISION = 4 

def prepKey (args): return "["+", ".join(args)+"]"

def fixedPrec(value):
   out = round(value*10**PRECISION)/10**PRECISION
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

def partition(V, CV1,CV2, EEV1,EEV2):
   CV = sorted(AA(sorted)(Delaunay(array(V)).vertices))
   BV1, BV2, BF1, BF2 = boundaryVertices( V, CV1,CV2, 'cuboid', EEV1,EEV2 )
   BV = BV1+BV2
   nE1 = len(EEV1)
   BF = BF1+[e+nE1 for e in BF2]
   return CV, BV1, BV2, BF1, BF2, BV, BF, nE1

""" Characteristic matrix transposition """
def invertRelation(V,CV):
   VC = [[] for k in range(len(V))]
   for k,cell in enumerate(CV):
      for v in cell:
         VC[v] += [k]
   return VC

""" Look for cells in Delaunay, with vertices in both operands """
def mixedCells(CV,CV1,CV2,n12):
   n0,n1 = 0, max(AA(max)(CV1))        # vertices in CV1 (extremes included)
   m0,m1 = n1+1-n12, max(AA(max)(CV2))    # vertices in CV2 (extremes included)
   return [list(cell) for cell in CV if any([ n0<=v<=n1 for v in cell]) 
      and any([ m0<=v<=m1 for v in cell])]

""" Look for cells in cells12, with vertices on boundaries """
def mixedCellsOnBoundaries(cells12,BV):
   cells12BV = [cell for cell in cells12
               if len(list(set(cell).intersection(BV))) != 0]
   return cells12BV

""" Build intersection tasks """
def cuttingTest(cuttingHyperplane,polytope,V):
   signs = [INNERPROD([cuttingHyperplane, V[v]+[1.]]) for v in polytope]
   signs = eval(vcode(signs))
   return any([value<-0.001 for value in signs]) and any([value>0.001 for value in signs])

def splittingTasks(V,pivots,BV,BC,VBC,CV,VC):
   tasks = []
   for pivotCell in pivots:
      cutVerts = [v for v in pivotCell if v in BV]
      for v in cutVerts:
         cutFacets = VBC[v]
         cells2cut = VC[v]
         for face,cell in CART([cutFacets,cells2cut]):
            polytope = CV[cell]
            points = [V[w] for w in BC[face]]
            dim = len(points[0])
            theMat = Matrix( [(dim+1)*[1.]] + [p+[1.] for p in points] )
            cuttingHyperplane = [(-1)**(col)*theMat.minor(0,col).determinant() 
                           for col in range(dim+1)]
            if cuttingTest(cuttingHyperplane,polytope,V):
               tasks += [[face,cell,cuttingHyperplane]]
   tasks = AA(eval)(set(AA(str)(tasks)))
   tasks = TrivialIntersection(tasks,V,BC,CV)
   return tasks

""" Trivial intersection filtering """
def TrivialIntersection(tasks,V,EEV,CV):
   out = []
   for face,cell,affineHull in tasks:
      faceVerts, cellVerts = EEV[face], CV[cell]
      v0 = list(set(faceVerts).intersection(cellVerts))[0] # v0 = common vertex
      transformMat = mat([VECTDIFF([V[v],V[v0]]) for v in cellVerts if v != v0]).T.I
      vects = (transformMat * (mat([VECTDIFF([V[v],V[v0]]) for v in faceVerts 
               if v != v0]).T)).T.tolist()
      if any([all([x>0 for x in list(vect)]) for vect in vects]): 
         out += [[face,cell,affineHull]]
   return out

""" Cell splitting in two cells """
def cellSplitting(face,cell,covector,V,EEV,CV):

   dim = len(V[0])
   subspace = (T(range(1,dim+1))(dim*[-50])(CUBOID(dim*[100])))
   normal = covector[:-1]
   if len(normal) == 2:  # 2D complex
      rotatedSubspace = R([1,2])(ATAN2(normal)-PI/2)(T(2)(-50)(subspace))
   elif len(normal) == 3:  # 3D complex
      rotatedSubspace = R()()(subspace)
   else: print "rotation error"
   t = V[EEV[face][0]]
   rototranslSubspace = T(range(1,dim+1))(t)(rotatedSubspace)
   cellHpc = MKPOL([V,[[v+1 for v in CV[cell]]],None])
   
   # cell1 = INTERSECTION([cellHpc,rototranslSubspace])
   tolerance=0.0001
   use_octree=False
   cell1 = Plasm.boolop(BOOL_CODE_AND, 
      [cellHpc,rototranslSubspace],tolerance,plasm_config.maxnumtry(),use_octree)
   verts,cells,pols = UKPOL(cell1)
   cell1 = AA(vcode)(verts)

   # cell2 = DIFFERENCE([cellHpc,rototranslSubspace]) 
   cell2 = Plasm.boolop(BOOL_CODE_DIFF, 
      [cellHpc,rototranslSubspace],tolerance,plasm_config.maxnumtry(),use_octree)
   verts,cells,pols = UKPOL(cell2)
   cell2 = AA(vcode)(verts)

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

""" Updating the split cell """
def splitCellUpdate(cell,vcell1,vcell2,CV):
   newVerts = list(set(vcell1).difference(CV[cell]))
   return newVerts

""" Updating the vertex set of split cells """
def splitCellsCreateVertices(vertdict,dict_fc,dict_cf,V,BC,CV,VC):
   DEBUG = False
   nverts = len(V); cellPairs = []; twoCellIndices = []; cuttingFaces = []
   while any([tasks != [] for face,tasks in dict_fc.items()]) : 
      for face,tasks in dict_fc.items():
         for task in tasks:
            cell,covector = task
            if cuttingTest(covector,CV[cell],V):
               cell1,cell2 = cellSplitting(face,cell,covector,V,BC,CV)
               if cell1 == [] or cell2 == []:
                  print "\nface,cell,covector =",face,cell,covector
                  print "cell1,cell2 =",cell1,cell2
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
                  V,CV, dict_cf, dict_fc,twoCells = splittingControl(face,cell,covector,vcell1,vcell2, 
                                          dict_fc,dict_cf,V,BC,CV,VC)
                  for adjCell in adjCells:
                     if cuttingTest(covector,CV[adjCell],V) and not ((face,covector) in dict_cf[adjCell]):
                        dict_fc[face] += [(adjCell,covector)] 
                        dict_cf[adjCell] += [(face,covector)] 
                  cellPairs += [[vcell1, vcell2]]
                  twoCellIndices += [twoCells]
                  cuttingFaces += [face]
               if DEBUG: showSplitting(V,cellPairs,BC,CV)
            else:
               dict_fc[face].remove((cell,covector))   # remove the split cell
               dict_cf[cell].remove((face,covector))   # remove the splitting face
   return cellPairs,twoCellIndices,cuttingFaces

""" Managing the splitting dictionaries """
def splittingControl(face,cell,covector,vcell1,vcell2,dict_fc,dict_cf,V,BC,CV,VC):

   print "vcell1,vcell2 =",vcell1,vcell2
   # only one facet covector crossing the cell
   cellVerts = CV[cell]
   CV[cell] = vcell1
   CV += [vcell2]
   twoCells = [cell,len(CV)-1]
   print "covector =",covector
   dict_fc[face].remove((cell,covector))   # remove the split cell
   dict_cf[cell].remove((face,covector))   # remove the splitting face
         
   # more than one facet covectors crossing the cell
   alist1,alist2 = list(),list()
   for aface,covector in dict_cf[cell]:
   
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
   
   dict_cf[cell] = alist1  
   dict_cf[len(CV)-1] = alist2
   return V,CV, dict_cf, dict_fc,twoCells

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

""" Splitting of Common Delaunay Complex """
def booleanBulk(V,n12,EEV,CV,VC,BF,CV1,CV2,EEV1,EEV2,BV,BV1,BV2,VEE1,VEE2):
   VE = [VEE1[v]+VEE2[v] for v in range(len(V))]
   cells12 = mixedCells(CV,CV1,CV2,n12)
   pivots = mixedCellsOnBoundaries(cells12,BV1,BV2)
   tasks = splittingTasks(V,pivots,BV,BF,VC,CV,EEV,VE)
      
   dict_fc,dict_cf = initTasks(tasks)
   vertdict = defaultdict(list)
   for k,v in enumerate(V): vertdict[vcode(v)] += [k]
   cellPairs,twoCellIndices = splitCellsCreateVertices(vertdict,dict_fc,dict_cf,V,EEV,CV,VC,BF)
   return cellPairs,twoCellIndices

""" Show the process of CDC splitting """
def showSplitting(V,cellPairs,BC,CV):
   VV = AA(LIST)(range(len(V)))
   boundaries = COLOR(RED)(SKEL_1(STRUCT(MKPOLS((V,BC)))))
   submodel = COLOR(CYAN)(STRUCT([ SKEL_1(STRUCT(MKPOLS((V,CV)))), boundaries ]))
   if cellPairs != []:
      cells1,cells2 = TRANS(cellPairs)
      out = [COLOR(WHITE)(MKPOL([V,[[v+1 for v in cell] for cell in cells1],None])), 
            COLOR(MAGENTA)(MKPOL([V,[[v+1 for v in cell] for cell in cells2],None]))]
      VIEW(STRUCT([ STRUCT(out),larModelNumbering(V,[VV,BC,CV],submodel,2) ]))
   else:
      VIEW(STRUCT([ larModelNumbering(V,[VV,BC,CV],submodel,2) ]))

""" Computation of bits of split cells """
def splitCellsBits(cuttingFaces,cellPairs,twoCellIndices,CV1,CV2,n12,BC):
   n0,n1 = 0, max(AA(max)(CV1))        # vertices in CV1 (extremes included)
   m0,m1 = n1+1-n12, max(AA(max)(CV2))    # vertices in CV2 (extremes included)
   print "\nn0,n1 =",n0,n1,"m0,m1 =",m0,m1
   for k,(v1,v2) in enumerate(BC):
      if v1>n1 or v2>n1: break
   boundarySpan1 = [0,k-1]
   boundarySpan2 = [k,len(BC)-1]
   print "boundarySpan1,boundarySpan2 =",boundarySpan1,boundarySpan2
   for face,cells,indices in zip(cuttingFaces,cellPairs,twoCellIndices):
      print "\ncells =", cells, "cell indices =", indices, "cutting face =",face
      cell1,cell2 = cells  # sets of vertex indices in V
      c1,c2 = indices  # cell indices in CV  (d-cells of CDC)
      v1s = list(set(cell1).difference(cell2))
      v2s = list(set(cell2).difference(cell1))
      faceVerts = BC[face]
      print "v1s,v2s =",v1s,v2s,"faceVerts =",faceVerts,
      if all([n0<=v<=n1 for v in v1s]) and all([m0<=v<=m1 for v in v2s]): print "bits = 1 0"
      elif all([n0<=v<=n1 for v in v2s]) and all([m0<=v<=m1 for v in v1s]): print "bits = 0 1"
      else: print "error"


