""" Module for Boolean ops with LAR """
DEBUG = True
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
global PRECISION
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


""" Characteristic matrix transposition """
def invertRelation(V,CV):
   VC = [[] for k in range(len(V))]
   for k,cell in enumerate(CV):
      for v in cell:
         VC[v] += [k]
   return VC





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


""" Updating the vertex set of split cells """
def testingSubspace(V,covector):
   def testingSubspace0(vcell):
      inout = SIGN(sum([INNERPROD([V[v]+[1.],covector]) for v in vcell]))
      return inout
   return testingSubspace0
   
def cuttingTest(cuttingHyperplane,polytope,V):
   signs = [INNERPROD([cuttingHyperplane, V[v]+[1.]]) for v in polytope]
   signs = eval(vcode(signs))
   return any([value<-0.001 for value in signs]) and any([value>0.001 for value in signs])

def splitCellsCreateVertices(vertdict,dict_fc,dict_cf,V,BC,CV,VC,lenBC1):
   CVbits = [[-1,-1] for k in range(len(CV))] 
   nverts = len(V); cellPairs = []; twoCellIndices = []; 
   while any([tasks != [] for face,tasks in dict_fc.items()]) : 
      for face,tasks in dict_fc.items():
         for task in tasks:
            cell,covector = task
            vcell = CV[cell]

            if cuttingTest(covector,vcell,V):
               cell1,cell2 = cellSplitting(face,cell,covector,V,BC,CV)
               if cell1 == [] or cell2 == []:
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
                                             
                  V,CV,CVbits, dict_cf, dict_fc,twoCells = splittingControl(
                     face,cell,covector,vcell,vcell1,vcell2, dict_fc,dict_cf,V,BC,CV,VC,CVbits,lenBC1)
                  if twoCells[0] != twoCells[1]:

                     for adjCell in adjCells:
                        dict_fc[face] += [(adjCell,covector)] 
                        dict_cf[adjCell] += [(face,covector)] 
                        cellPairs += [[vcell1, vcell2]]
                        twoCellIndices += [[twoCells]]
                                    
               DEBUG = False
               if DEBUG: showSplitting(V,cellPairs,BC,CV)
            else:
               dict_fc[face].remove((cell,covector))   # remove the split cell
               dict_cf[cell].remove((face,covector))   # remove the splitting face
   return CVbits,cellPairs,twoCellIndices

""" Managing the splitting dictionaries """
def splittingControl(face,cell,covector,vcell,vcell1,vcell2,dict_fc,dict_cf,V,BC,CV,VC,CVbits,lenBC1):

   boundaryFacet = BC[face]
   translVector = V[boundaryFacet[0]]
   tcovector = [cv+tv*covector[-1] for (cv,tv) in zip(covector[:-1],translVector) ]+[0.0]

   c1,c2 = cell,cell
   if not haltingSplitTest(cell,vcell,vcell1,vcell2,boundaryFacet,translVector,tcovector,V) :

      # only one facet covector crossing the cell
      cellVerts = CV[cell]
      CV[cell] = vcell1
      CV += [vcell2]
      CVbits += [copy(CVbits[cell])]
      c1,c2 = cell,len(CV)-1
   
      firstCell,secondCell = AA(testingSubspace(V,covector))([vcell1,vcell2])
      if face < lenBC1 and firstCell==-1:       # face in boundary(op1)
         CVbits[c1][0] = 0
         CVbits[c2][0] = 1
      elif face >= lenBC1 and firstCell==-1:    # face in boundary(op2)
         CVbits[c1][1] = 0 
         CVbits[c2][1] = 1
      else: print "error splitting face,c1,c2 =",face,c1,c2
   
      #dict_fc[face].remove((cell,covector)) # remove the split cell
      #dict_cf[cell].remove((face,covector)) # remove the splitting face
            
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
      
   else:
      dict_fc[face].remove((cell,covector))  # remove the split cell
      dict_cf[cell].remove((face,covector))  # remove the splitting face   
      
   return V,CV,CVbits, dict_cf, dict_fc,[c1,c2]

""" Test for split halting along a boundary facet """
def haltingSplitTest(cell,vcell,vcell1,vcell2,boundaryFacet,translVector,tcovector,V):
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
   if verts == []: 
      print "\n****** cell =",cell
      return True
   else: return False

# cell1 = INTERSECTION([cellHpc,rototranslSubspace])
# tolerance=0.0001
# use_octree=False
# cell1 = Plasm.boolop(BOOL_CODE_AND, 
#  [cellHpc,rototranslSubspace],tolerance,plasm_config.maxnumtry(),use_octree)
# verts,cells,pols = UKPOL(cell1)
# cell1 = AA(vcode)(verts)
# if 

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
   CV = sorted(AA(sorted)(Delaunay(array(V)).vertices))
   vertdict = defaultdict(list)
   for k,v in enumerate(V): vertdict[vcode(v)] += [k]
   
   BC1 = signedCellularBoundaryCells(V1,basis1)
   BC2 = signedCellularBoundaryCells(V2,basis2)
   BC = sorted([[ vertdict[vcode(V1[v])][0] for v in cell] for cell in BC1] + [ 
         [ vertdict[vcode(V2[v])][0] for v in cell] for cell in BC2])
   BV = list(set(CAT([v for v in BC])))
   VV = AA(LIST)(range(len(V)))

   print "\n BC =",BC,'\n'


   if DEBUG: 
      """ Input and CDC visualisation """
      submodel1 = mkSignedEdges((V1,BC1))
      submodel2 = mkSignedEdges((V2,BC2))
      VIEW(STRUCT([submodel1,submodel2]))
      submodel = SKEL_1(STRUCT(MKPOLS((V,CV))))
      VIEW(larModelNumbering(V,[VV,BC,CV],submodel,4))
      submodel = STRUCT([SKEL_1(STRUCT(MKPOLS((V,CV)))), COLOR(RED)(STRUCT(MKPOLS((V,BC))))])
      VIEW(larModelNumbering(V,[VV,BC,CV],submodel,4))
      
   """ New implementation of splitting dictionaries """
   VC = invertRelation(V,CV)
   
   covectors = []
   for faceVerts in BC:
      points = [V[v] for v in faceVerts]
      dim = len(points[0])
      theMat = Matrix( [(dim+1)*[1.]] + [p+[1.] for p in points] )
      covector = [(-1)**(col)*theMat.minor(0,col).determinant() 
                     for col in range(dim+1)]
      covectors += [covector]
   
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
         else: print "error: found no face,cell,covector"
   
   """ Initialization of splitting dictionaries """
   tasks = []
   for face,covector in zip(range(len(BC)),covectors):
      tasks += [covectorCell(face,BC[face],covector,CV,VC)]
   
   dict_fc,dict_cf = initTasks(tasks)
   
   
   
   CVbits,cellPairs,twoCellIndices = splitCellsCreateVertices( 
      vertdict,dict_fc,dict_cf,V,BC,CV,VC,len(BC1))
   showSplitting(V,cellPairs,BC,CV)
   
   print "\n"
   for k in range(len(CV)):  print "k,CVbits[k],CV[k] =",k,CVbits[k],CV[k]
   
   for cell in range(len(CV)):
      if CVbits[cell][0] == 1:
         CVbits = booleanChainTraverse(0,cell,V,CV,CVbits,1)      
      if CVbits[cell][0] == 0:
         CVbits = booleanChainTraverse(0,cell,V,CV,CVbits,0)
      if CVbits[cell][1] == 1:
         CVbits = booleanChainTraverse(1,cell,V,CV,CVbits,1)
      if CVbits[cell][1] == 0:
         CVbits = booleanChainTraverse(1,cell,V,CV,CVbits,0)
   
   print "\n"
   for k in range(len(CV)):  print "k,CVbits[k],CV[k] =",k,CVbits[k],CV[k]
   
   chain1,chain2 = TRANS(CVbits)
   print "\ndict_cf",dict_cf
   print "\ndict_fc",dict_fc,"\n"
   return V,CV,chain1,chain2,CVbits

