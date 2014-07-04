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
         for facet,cell2cut in CART([cutFacets,cells2cut]):
            polytope = CV[cell2cut]
            points = [V[w] for w in EEV[facet]]
            dim = len(points[0])
            theMat = Matrix( [(dim+1)*[1.]] + [p+[1.] for p in points] )
            cuttingHyperplane = [(-1)**(col)*theMat.minor(0,col).determinant() 
                           for col in range(dim+1)]
            if cuttingTest(cuttingHyperplane,polytope,V):
               tasks += [[facet,cell2cut,cuttingHyperplane]]
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
      transformMat = mat([VECTDIFF([V[v],V[v0]]) for v in cellVerts if v != v0]).I
      vects = (transformMat * mat([VECTDIFF([V[v],V[v0]]) for v in faceVerts 
               if v != v0]).T).tolist()
      if any([all([x>0 for x in list(vect)]) for vect in vects]): out += [[face,cell,affineHull]]
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
   print "\ncell =",UKPOL(cellHpc)[0]
   # cell1 = INTERSECTION([cellHpc,rototranslSubspace])
   
   tolerance=0.0001
   use_octree=False
   cell1 = Plasm.boolop(BOOL_CODE_AND, 
      [cellHpc,rototranslSubspace],tolerance,plasm_config.maxnumtry(),use_octree)
   
   # cell2 = DIFFERENCE([cellHpc,rototranslSubspace])
   
   cell2 = Plasm.boolop(BOOL_CODE_DIFF, 
      [cellHpc,rototranslSubspace],tolerance,plasm_config.maxnumtry(),use_octree)

   return cell1,cell2

