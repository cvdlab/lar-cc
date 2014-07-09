""" Module for Boolean ops with LAR """
from matrix import *
from pyplasm import *
from scipy import *
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


def covering(model1,model2):
   V, CV1, CV2, n12 = vertexSieve(model1,model2)
   _,EEV1 = larFacets((V,CV1),dim=2,emptyCellNumber=1)
   _,EEV2 = larFacets((V,CV2),dim=2,emptyCellNumber=1)
   CV1 = CV1[:-1]
   CV2 = CV2[:-1]
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
def mixedCells(CV,n0,n1,m0,m1):
   return [list(cell) for cell in CV if any([ n0<=v<=n1 for v in cell]) 
      and any([ m0<=v<=m1 for v in cell])]

""" Look for cells in cells12, with vertices on boundaries """
def mixedCellsOnBoundaries(cells12,BV1,BV2):
   cells12BV1 = [cell for cell in cells12
               if len(list(set(cell).intersection(BV1))) != 0]
   cells12BV2 = [cell for cell in cells12
               if len(list(set(cell).intersection(BV2))) != 0]
   pivots = sorted(AA(sorted)(cells12BV1+cells12BV2))
   pivots = [cell for k,cell in enumerate(pivots[:-1]) if cell==pivots[k+1]]
   return pivots

""" Build intersection tasks """
def cuttingTest(cuttingHyperplane,polytope,V):
   signs = [INNERPROD([cuttingHyperplane, V[v]+[1.]]) for v in polytope]
   return any([value<0 for value in signs]) and any([value>0 for value in signs])

def splittingTasks(V,pivots,BV,BF,VC,CV,EEV,VE):
   tasks = []
   for pivotCell in pivots:
      cutVerts = [v for v in pivotCell if v in BV]
      for v in cutVerts:
         cutFacets = [e for e in VE[v] if e in BF]
         cells2cut = VC[v]
         for face,cell in CART([cutFacets,cells2cut]):
            polytope = CV[cell]
            points = [V[w] for w in EEV[face]]
            dim = len(points[0])
            theMat = Matrix( [(dim+1)*[1.]] + [p+[1.] for p in points] )
            cuttingHyperplane = [(-1)**(col)*theMat.minor(0,col).determinant() 
                           for col in range(dim+1)]
            if cuttingTest(cuttingHyperplane,polytope,V):
               tasks += [[face,cell,cuttingHyperplane]]
   tasks = AA(eval)(set(AA(str)(tasks)))
   tasks = TrivialIntersection(tasks,V,EEV,CV)
   print "\ntasks =",tasks
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
   print "\nV,CV =",V,CV
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
   print "cell,vcell1,vcell2 =",cell,vcell1,vcell2
   newVerts = list(set(vcell1).difference(CV[cell]))
   return newVerts

""" Updating the vertex set of split cells """
def splitCellsCreateVertices(dict_fc,dict_cf,V,EEV,CV,VC):
   out = []; nverts = len(V); cellPairs = []
   vertdict = defaultdict(list)
   for k,v in enumerate(V): vertdict[vcode(v)] += [k]
   while any([item[1]!=[] for item in dict_fc.items()]) : 
      for face,tasks in dict_fc.items():
         for task in tasks:
            cell,covector = task
            print "\nface,cell,covector =",face,cell,covector
            cell1,cell2 = cellSplitting(face,cell,covector,V,EEV,CV)
            print "\ncell1,cell2 =",cell1,cell2
            vcell1 = []
            for k in cell1:
               if vertdict[k]==[]: 
                  vertdict[k] += [nverts]
                  V += [eval(k)]
                  nverts += 1
               vcell1 += [vertdict[k]]
            vcell1 = CAT(vcell1)
            vcell2 = CAT([vertdict[k] for k in cell2])
            newVerts = splitCellUpdate(cell,vcell1,vcell2,CV)
            print "\nnewVerts =",newVerts
            V,CV, dict_cf, dict_fc = splittingControl(face,cell,vcell1,vcell2, 
                                    dict_fc,dict_cf,V,EEV,CV,VC)
            cellPairs += [[vcell1, vcell2]]
   print "\n***** dict_fc.items()",dict_fc.items()
   print "\n***** dict_cf.items()",dict_cf.items()
   return vertdict, cellPairs, nverts

""" Managing the splitting dictionaries """
def splittingControl(face,cell,vcell1,vcell2,dict_fc,dict_cf,V,EEV,CV,VC):

   # only one facet covector crossing the cell

   print "***** 1"
   cellVerts = CV[cell]
   print "cellVerts =",cellVerts
   CV[cell] = vcell1
   CV += [vcell2]
   covector = dict_cf[cell][0][1]
   dict_fc[face].remove((cell,covector))   # remove the split cell
   dict_cf[cell].remove((face,covector))   # remove the splitting face
   print "dict_fc =",dict_fc
   print "dict_cf =",dict_cf
         
   # more than one facet covectors crossing the cell
      
   print "***** 2"
   alist1,alist2 = list(),list()
   print "alist1,alist2 =",alist1,alist2
   
   for aface,covector in dict_cf[cell]:
   
      # for each facet crossing the cell
      # compute the intersection between the facet and the cell
      
      print "***** 3"
      print "aface,covector =",aface,covector
            
      faceVerts = EEV[aface]
      print "faceVerts =",faceVerts
      commonVerts = list(set(faceVerts).intersection(cellVerts))
      print "commonVerts =",commonVerts
      
      # and attribute the intersection to the split subcells
      
      if set(vcell1).intersection(commonVerts) != set():
         alist1.append((aface,covector))
         print "alist1 =",alist1
      else: dict_fc[aface].remove((cell,covector)) 
            
      if set(vcell2).intersection(commonVerts) != set():
         alist2.append((aface,covector))
         print "alist2 =",alist2
         dict_fc[aface] += [(len(CV)-1,covector)]
   
   print "***** 4"
   dict_cf[cell] = alist1  
   dict_cf[len(CV)-1] = alist2
   print "dict_cf =",dict_cf
   # for f in alist1: dict_fc[f] = cell
   # for f in alist2: dict_fc[f] = len(CV)-1

   return V,CV, dict_cf, dict_fc

