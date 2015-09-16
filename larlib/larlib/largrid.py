"""Module with functions for grid generation and Cartesian product"""
import collections
from larlib import *

def larSplit(dom):
    def larSplit1(n):
        # assert n > 0 and isinstance(n,int)
        item = float(dom)/n
        ints = range(n+1)
        items = [item]*(n+1)
        vertices = [[int*item] for (int,item) in zip(ints,items)]
        return vertices
    return larSplit1

def grid0(n):
    cells = AA(LIST)(range(n+1))
    return cells

def grid1(n):
    ints = range(n+1)
    cells = TRANS([ints[:-1],ints[1:]])
    return cells

def larGrid(n):
    def larGrid1(d):
        if d==0: return grid0(n)
        elif d==1: return grid1(n)
    return larGrid1

def larVertProd(vertLists):
    return AA(CAT)(CART(vertLists))

def index2addr (shape):
    n = len(shape)
    shape = shape[1:]+[1]
    weights = [PROD(shape[k:]) for k in range(n)]
    def index2addr0 (multindex):
        return INNERPROD([multindex, weights])
    return index2addr0

def larCellProd(cellLists):
    shapes = [len(item) for item in cellLists]
    indices = CART([range(shape) for shape in shapes])
    jointCells = [CART([cells[k] for k,cells in zip(index,cellLists)])
                  for index in indices]
    convert = index2addr([ shape+1 if (len(cellLists[k][0]) > 1) else shape
                             for k,shape in enumerate(shapes) ])
    return [AA(convert)(cell) for cell in jointCells]

def binaryRange(n):
    return [('{0:0'+str(n)+'b}').format(k) for k in range(2**n)]

def filterByOrder(n):
    terms = [AA(int)(list(term)) for term in binaryRange(n)]
    return [[term for term in terms if sum(term) == k] for k in range(n+1)]

def larGridSkeleton(shape):
    n = len(shape)
    def larGridSkeleton0(d):
        components = filterByOrder(n)[d]
        componentCellLists = [AA(APPLY)(zip( AA(larGrid)(shape),(component) ))
                              for component in components]
        return CAT([ larCellProd(cellLists)  for cellLists in componentCellLists ])
    return larGridSkeleton0

def larImageVerts(shape):
   def vertexDomain(n): 
      return [[k] for k in range(n)]
   vertLists = [vertexDomain(k+1) for k in shape]
   vertGrid = larVertProd(vertLists)
   return vertGrid

def larCuboids(shape, full=False):
   vertGrid = larImageVerts(shape)
   gridMap = larGridSkeleton(shape)
   if not full: 
      cells = gridMap(len(shape))
   else:
      skeletonIds = range(len(shape)+1)
      cells = [ gridMap(id) for id in skeletonIds ]
   return vertGrid, cells

def gridSkeletons(shape):
   gridMap = larGridSkeleton(shape)
   skeletonIds = range(len(shape)+1)
   skeletons = [ gridMap(id) for id in skeletonIds ]
   return skeletons
   
if __name__=="__main__":
   print "\ngridSkeletons([3]) =\n", gridSkeletons([3])
   print "\ngridSkeletons([3,2]) =\n", gridSkeletons([3,2])
   print "\ngridSkeletons([3,2,1]) =\n", gridSkeletons([3,2,1])

def gridBoundaryMatrices(shape):
   skeletons = gridSkeletons(shape)
   boundaryMatrices = [boundary(skeletons[k+1],faces) 
                   for k,faces in enumerate(skeletons[:-1])]
   return boundaryMatrices
   
if __name__=="__main__":
   for k in range(1):
      print "\ngridBoundaryMatrices([3]) =\n", \
            csr2DenseMatrix(gridBoundaryMatrices([3])[k])
   for k in range(2):
      print "\ngridBoundaryMatrices([3,2]) =\n", \
            csr2DenseMatrix(gridBoundaryMatrices([3,2])[k])
   for k in range(3):
      print "\ngridBoundaryMatrices([3,2,1]) =\n", \
            csr2DenseMatrix(gridBoundaryMatrices([3,2,1])[k])

def larModelProduct(twoModels):
    (V, cells1), (W, cells2) = twoModels
    vertices = collections.OrderedDict(); k = 0
    for v in V:
        for w in W:
            id = tuple(v+w)
            if not vertices.has_key(id):
                vertices[id] = k
                k += 1   
    cells = [ [vertices[tuple(V[v] + W[w])] for v in c1 for w in c2]
             for c1 in cells1 for c2 in cells2]  
    model = [list(v) for v in vertices.keys()], cells
    return model

""" Simplicial face stack computation """
def larSimplicialStack(simplices):
   dim = len(simplices[0])-1
   faceStack = [simplices]
   for k in range(dim):
      faces = larSimplexFacets(faceStack[-1])
      faceStack.append(faces)
   return REVERSE(faceStack)

""" Extraction of facets from cuboidal complexes """
def larCuboidsFacets((V,cells)):
   dim = len(V[0])
   n = int(2**(dim-1))
   facets = []
   for cell in cells:
      coords = [AR([V[v],v]) for v in cell] # decorate coords with vertex index
      doubleFacets = [sorted(coords,key=(lambda x: x[k])) for k in range(dim)]
      facets += AA(AA(LAST))(CAT([[pair[:n],pair[n:]] for pair in doubleFacets]))
   facets = AA(eval)(set(AA(str)(facets))) # remove duplicates
   return V,sorted(facets)

if __name__ == "__main__":
   VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(larCuboidsFacets(larCuboids([3,3,3])))))

if __name__=="__main__":
   def mergeSkeletons(larSkeletons): return larSkeletons[0],CAT(larSkeletons[1])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(mergeSkeletons(larCuboids([3],True)))))
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(mergeSkeletons(larCuboids([3,2],True)))))
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(mergeSkeletons(larCuboids([3,2,1],True)))))
   
   if __name__ == "__main__":
       geom_0,topol_0 = [[0.],[1.],[2.],[3.],[4.]],[[0,1],[1,2],[2,3],[3,4]]
       geom_1,topol_1 = [[0.],[1.],[2.]], [[0,1],[1,2]]
       mod_0 = (geom_0,topol_0)
       mod_1 = (geom_1,topol_1)
       squares = larModelProduct([mod_0,mod_1])
       VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(squares)))
       cubes = larModelProduct([squares,mod_0])
       VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(cubes)))
   
