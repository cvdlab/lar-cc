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

from splitcell import *
DEBUG = False
""" TODO: use package Decimal (http://docs.python.org/2/library/decimal.html) """
global PRECISION
PRECISION = 4.

def verySmall(number): return abs(number) < 10**-(PRECISION)

def prepKey (args): return "["+", ".join(args)+"]"

def fixedPrec(value):
   out = round(value*10**(PRECISION))/10**(PRECISION)
   if out == -0.0: out = 0.0
   return str(out)
   
def vcode (vect): 
   """
   To generate a string representation of a number array.
   Used to generate the vertex keys in PointSet dictionary, and other similar operations.
   """
   return prepKey(AA(fixedPrec)(vect))

""" Merge two dictionaries with keys the point locations """
def mergeVertices(model1, model2):

   (V1,CV1),(V2,CV2) = model1, model2

   n = len(V1); m = len(V2)
   def shift(CV, n): 
      return [[v+n for v in cell] for cell in CV]
   CV2 = shift(CV2,n)

   vdict1 = defaultdict(list)
   for k,v in enumerate(V1): vdict1[vcode(v)].append(k) 
   vdict2 = defaultdict(list)
   for k,v in enumerate(V2): vdict2[vcode(v)].append(k+n) 
   vertDict = defaultdict(list)
   for point in vdict1.keys(): vertDict[point] += vdict1[point]
   for point in vdict2.keys(): vertDict[point] += vdict2[point]

   case1, case12, case2 = [],[],[]
   for item in vertDict.items():
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

   return V,CV1,CV2, n1+n2,n2,n2+n3

""" Make Common Delaunay Complex """
def makeCDC(arg1,arg2):

   (V1,basis1), (V2,basis2) = arg1,arg2
   (facets1,cells1),(facets2,cells2) = basis1[-2:],basis2[-2:]
   model1, model2 = (V1,cells1),(V2,cells2)

   V, _,_, n1,n12,n2 = mergeVertices(model1, model2)
   n = len(V)
   assert n == n1 - n12 + n2
   
   CV = sorted(AA(sorted)(Delaunay(array(V)).simplices))
   
   vertDict = defaultdict(list)
   for k,v in enumerate(V): vertDict[vcode(v)] += [k]
   
   BC1 = signedCellularBoundaryCells(V1,basis1)
   BC2 = signedCellularBoundaryCells(V2,basis2)
   BC = [[ vertDict[vcode(V1[v])][0] for v in cell] for cell in BC1] + [ 
         [ vertDict[vcode(V2[v])][0] for v in cell] for cell in BC2] #+ qhullBoundary(V)
   
   return V,CV,vertDict,n1,n12,n2,BC,len(BC1),len(BC2)

""" Cell-facet intersection test """
def cellFacetIntersecting(boundaryFacet,cell,covector,V,CV):
   points = [V[v] for v in CV[cell]]
   vcell1,newFacet,vcell2 = SPLITCELL(covector,points,tolerance=1e-4,ntry=4)
   boundaryFacet = [V[v] for v in boundaryFacet]
   translVector = boundaryFacet[0]
   
   # translation 
   newFacet = [ VECTDIFF([v,translVector]) for v in newFacet ]
   boundaryFacet = [ VECTDIFF([v,translVector]) for v in boundaryFacet ]
   
   # linear transformation: boundaryFacet -> standard (d-1)-simplex
   d = len(V[0])
   transformMat = mat( boundaryFacet[1:d] + [covector[1:]] ).T.I
   
   # transformation in the subspace x_d = 0
   newFacet = (transformMat * (mat(newFacet).T)).T.tolist()
   boundaryFacet = (transformMat * (mat(boundaryFacet).T)).T.tolist()
   
   # projection in E^{d-1} space and Boolean test
   newFacet = MKPOL([ AA(lambda v: v[:-1])(newFacet), [range(1,len(newFacet)+1)], None ])
   boundaryFacet = MKPOL([ AA(lambda v: v[:-1])(boundaryFacet), [range(1,len(boundaryFacet)+1)], None ])
   verts,cells,pols = UKPOL(INTERSECTION([newFacet,boundaryFacet]))
   
   if verts == []: return False
   else: return True


""" Splitting tests """
def testingSubspace(V,covector):
   def testingSubspace0(vcell):
      inout = SIGN(sum([INNERPROD([[1.]+V[v],covector]) for v in vcell]))
      return inout
   return testingSubspace0
   
def cuttingTest(covector,polytope,V):
   signs = [INNERPROD([covector, [1.]+V[v]]) for v in polytope]
   signs = eval(vcode(signs))
   return any([value<-0.001 for value in signs]) and any([value>0.001 for value in signs])
   
def tangentTest(covector,facet,adjCell,V):
   common = list(set(facet).intersection(adjCell))
   signs = [INNERPROD([covector, [1.]+V[v]]) for v in common]
   count = 0
   for value in signs:
      if -0.0001<value<0.0001: count +=1
   if count >= len(V[0]): 
      return True
   else: 
      return False   


""" Elementary splitting test """
def dividenda(V,CV, cell,facet,covector,unchosen):
   out = []
   adjCells = adjacencyQuery(V,CV)(cell)
   for adjCell in set(adjCells).difference(unchosen):
      if (cuttingTest(covector,CV[adjCell],V) and \
         cellFacetIntersecting(facet,adjCell,covector,V,CV)) or \
         tangentTest(covector,facet,CV[adjCell],V): out += [adjCell]
   return out

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

""" Computation of boundary facets covering with CDC cells """
def boundaryCover(V,CV,BC,VC):
   cellsToSplit = list()
   boundaryCellCovering = []
   for k,facet in enumerate(BC):
      covector = COVECTOR([V[v] for v in facet])
      seedsOnFacet = VC[facet[0]]
      cellsToSplit = [dividenda(V,CV, cell,facet,covector,[]) for cell in seedsOnFacet ]
      cellsToSplit = set(CAT(cellsToSplit))
      while True:
         newCells = [dividenda(V,CV, cell,facet,covector,cellsToSplit) for cell in cellsToSplit ]
         if newCells != []: newCells = CAT(newCells)
         covering = cellsToSplit.union(newCells)
         if covering == cellsToSplit: 
            break
         cellsToSplit = covering
      boundaryCellCovering += [list(covering)]
   return boundaryCellCovering

""" CDC cell splitting with one or more cutting facets """
def fragment(cell,cellCuts,V,CV,BC):
   vcell = CV[cell]
   cellFragments = [[V[v] for v in vcell]]
   
   for f in cellCuts[cell]:
      facet = BC[f]
      plane = COVECTOR([V[v] for v in facet])
      for k,fragment in enumerate(cellFragments):
      
         #if not tangentTest(plane,facet,fragment,V):
         [below,equal,above] = SPLITCELL(plane,fragment,tolerance=1e-4,ntry=4)
         if below != above:
            cellFragments[k] = below
            cellFragments += [above]
      facets = facetsOnCuts(cellFragments,cellCuts,V,BC)
   return cellFragments

""" SCDC splitting with every boundary facet """
def makeSCDC(V,CV,BC,nbc1,nbc2):
   index,defaultValue = -1,-1
   VC = invertRelation(CV)
   CW,BCfrags = [],[]
   Wdict = dict()
   BCellcovering = boundaryCover(V,CV,BC,VC)

   cellCuts = invertRelation(BCellcovering)
   for k in range(len(CV) - len(cellCuts)): cellCuts += [[]]
   
   def verySmall(number): return abs(number) < 10**-5.5
   
   for k,frags in enumerate(cellCuts):
      if cellCuts[k] == []:
         cell = []
         for v in CV[k]:
            key = vcode(V[v])
            if Wdict.get(key,defaultValue) == defaultValue:
               index += 1
               Wdict[key] = index
               cell += [index]
            else: 
               cell += [Wdict[key]]
         CW += [cell]
      else:
         cellFragments = fragment(k,cellCuts,V,CV,BC)
         for cellFragment in cellFragments:
            cellFrag = []
            for v in cellFragment:
               key = vcode(v)
               if Wdict.get(key,defaultValue) == defaultValue:
                  index += 1
                  Wdict[key] = index
                  cellFrag += [index]
               else: 
                  cellFrag += [Wdict[key]]
            CW += [cellFrag]  
            
            BCfrags += [ (h, [Wdict[vcode(w)] for w in cellFragment if verySmall( 
                        PROD([ COVECTOR( [V[v] for v in BC[h]] ), [1.]+w ])) ] )
                      for h in cellCuts[k]]  
   
   BCW = [ [ Wdict[vcode(V[v])] for v in cell ] for cell in BC]
         
         
   W = sorted(zip( Wdict.values(), Wdict.keys() ))
   W = AA(eval)(TRANS(W)[1])
   dim = len(W[0])
   
   boundary1,boundary2 = defaultdict(list),defaultdict(list)                   
   for h,frags in BCfrags:
      if h < nbc1: boundary1[h] += [frags]
      else: boundary2[h] += [frags] 
   boundarylist1,boundarylist2 = [],[]
   for h,facets in boundary1.items():
      boundarylist1 += [(h, AA(eval)(set([str(sorted(f)) for f in facets if len(set(f)) >= dim])) )]
   for h,facets in boundary2.items():
      boundarylist2 += [(h, AA(eval)(set([str(sorted(f)) for f in facets if len(set(f)) >= dim])) )]
   boundary1,boundary2 = dict(boundarylist1),dict(boundarylist2)

   return W,CW,VC,BCellcovering,cellCuts,boundary1,boundary2,BCW

""" Characteristic matrix transposition """
def invertRelation(CV):
   def myMax(List):
      if List==[]: return -1
      else: return max(List)
   columnNumber = max(AA(myMax)(CV))+1
   VC = [[] for k in range(columnNumber)]
   for k,cell in enumerate(CV):
      for v in cell:
         VC[v] += [k]
   return VC

""" Computation of embedded boundary cells """
def facetsOnCuts(cellFragments,cellCuts,V,BC):


   pass
   return #facets

""" Coboundary operator on the convex decomposition of common space """
from scipy.spatial import ConvexHull
def qhullBoundary(V):
   points = array(V)
   hull = ConvexHull(points)
   out = hull.simplices.tolist()
   return sorted(out)

""" Extracting a $(d-1)$-basis of SCDC """
def larConvexFacets (V,CV):
   dim = len(V[0])
   model = V,CV
   V,FV = larFacets(model,dim)
   FV = AA(eval)(list(set(AA(str)(FV + convexBoundary(V,CV)))))
   return sorted(AA(sorted)(FV))

""" Computation of boundary operator of a convex LAR model"""
def convexBoundary(W,CW):
   points = array(W)
   hull = ConvexHull(points,qhull_options="Qc")
   coplanarVerts = hull.coplanar.tolist()
   if coplanarVerts != []:  coplanarVerts = CAT(coplanarVerts)
   BWchain = set( CAT(qhullBoundary(W)) + coplanarVerts )
   dim = len(W[0])
   bfacets = [list(BWchain.intersection(cell)) 
               for cell in CW if len(BWchain.intersection(cell)) >= dim]
   return bfacets

""" Writing labelling seeds on SCDC """
def cellTagging(boundaryDict,boundaryMat,CW,FW,W,BC,CWbits,arg):
   dim = len(W[0])
   for face in boundaryDict:
      for facet in boundaryDict[face]:
         cofaces = list(boundaryMat[facet].tocoo().col)
         if len(cofaces) == 1: 
            CWbits[cofaces[0]][arg] = 1
         elif len(cofaces) == 2:
            v0 = list(set(CW[cofaces[0]]).difference(FW[facet]))[0]
            v1 = list(set(CW[cofaces[1]]).difference(FW[facet]))[0]
            # take d affinely independent vertices in face (TODO: use pivotSimplices() 
            simplex0 = BC[face][:dim] + [v0]
            simplex1 = BC[face][:dim] + [v1]
            sign0 = sign(det([W[v]+[1] for v in simplex0]))
            sign1 = sign(det([W[v]+[1] for v in simplex1]))
            
            if sign0 == 1: CWbits[cofaces[0]][arg] = 1
            elif sign0 == -1: CWbits[cofaces[0]][arg] = 0
            if sign1 == 1: CWbits[cofaces[1]][arg] = 1
            elif sign1 == -1: CWbits[cofaces[1]][arg] = 0
         else: 
            print "error: too many cofaces of boundary facets"
   return CWbits

""" Recursive diffusion of labels on SCDC """
def booleanChainTraverse(h,cell,V,CV,CWbits,value):
   adjCells = adjacencyQuery(V,CV)(cell)
   for adjCell in adjCells: 
      if CWbits[adjCell][h] == -1:
         CWbits[adjCell][h] = value
         CWbits = booleanChainTraverse(h,adjCell,V,CV,CWbits,value)
   return CWbits

""" Mapping from hyperplanes to lists of facets """
def facet2covectors(W,FW):
   return [COVECTOR([W[v] for v in facet]) for facet in FW]

def boundaries(boundary1,boundary2):
   return set(CAT(boundary1.values() + boundary2.values()))

""" Mapping from hyperplanes to lists of facets """
def facet2covectors(W,FW):
   return [COVECTOR([W[v] for v in facet]) for facet in FW]

def boundaries(boundary1,boundary2):
   return set(CAT(boundary1.values() + boundary2.values()))

""" Building the boundary complex of the current chain """
def chain2complex(W,CW,FW,chain,coBoundaryMat,boundaryMat,constraints):
   chainCoords = csc_matrix((len(CW), 1))
   for cell in chain: chainCoords[cell] = 1
   out = boundaryCells = set((boundaryMat * chainCoords
                     ).tocoo().row).difference(constraints)
   return out

""" Sticking cells together """
def testAttachment(usedCells,facet,chain, 
         W,FW,coBoundaryMat,envelope,covectors):
   pass

def protrudeChain (W,CW,FW,chain,coBoundaryMat,boundaryMat,
                  envelope,covectors,usedCells):
   for facet in envelope:
      done,usedCells = testAttachment(usedCells,facet,chain, 
         W,FW,coBoundaryMat,envelope,covectors)
      if done: 
         usedCells[facet] = True
         verts = list(set([FW[cell] for cell in envelope+[facet]]))
         break
   return done,verts

""" Gathering and writing a polytopal complex """
def gatherPolytopes(W,CW,FW,coBoundaryMat,boundaryMat,covectors,
                  bounds1,bounds2):
   usedCells = [False for cell in CW]
   covectors = facet2covectors(W,FW)
   constraints = boundaries(bounds1,bounds2)
   Xdict,index = dict(),0
   while not all(usedCells):
      for cell in CW:
         if not usedCells[cell]:
            chain = [cell]
            usedCells[cell] = True
            done,verts = protrudeChain (W,CW,FW,chain,coBoundaryMat,boundaryMat,
                  envelope,covectors,usedCells)
            if done:
               for v in verts:
                  Xdict[vcode(W[v])] = index
                  index += 1
               CX += [chain]
      X = AA(eval)(TRANS(sorted(Xdict.values(),Xdict.keys()))[1])
   return X,CX


