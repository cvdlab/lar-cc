from larcc import *
from lar2psm import *
from scipy import *

#------------------------------------------------------------------
# extrusion of simplicial complexes
# combinatorial algorithm

def VERTEXTRUDE((V,coords)):
    """
        Utility function to generate the output model vertices in a multiple extrusion of a LAR model.
        V is a list of d-vertices (each given as a list of d coordinates).
        coords is a list of absolute translation parameters to be applied to V in order
        to generate the output vertices.
        
        Return a new list of (d+1)-vertices.
        """
    return CAT(AA(COMP([AA(AR),DISTR]))(DISTL([V,coords])))

def cumsum(iterable):
    # cumulative addition: list(cumsum(range(4))) => [0, 1, 3, 6]
    iterable = iter(iterable)
    s = iterable.next()
    yield s
    for c in iterable:
        s = s + c
        yield s

def larExtrude(model,pattern):
    V,FV = model
    d = len(FV[0])
    offset = len(V)
    m = len(pattern)
    outcells = []
    for cell in FV:
        # create the indices of vertices in the cell "tube"
        tube = [v + k*offset for k in range(m+1) for v in cell]
        # take groups of d+1 elements, via shifting by one
        rangelimit = len(tube)-d
        cellTube = [tube[k:k+d+1] for k in range(rangelimit)]
        outcells += [scipy.reshape(cellTube,newshape=(m,d,d+1)).tolist()]
    outcells = AA(CAT)(TRANS(outcells))
    outcells = [group for k,group in enumerate(outcells) if pattern[k]>0 ]
    coords = list(cumsum([0]+(AA(ABS)(pattern))))
    outVerts = VERTEXTRUDE((V,coords))
    newModel = outVerts, CAT(outcells)
    return newModel

if __name__ == "__main__":
    V = [[0,0],[1,0],[2,0],[0,1],[1,1],[2,1],[0,2],[1,2],[2,2]]
    FV = [[0,1,3],[1,2,4],[2,4,5],[3,4,6],[4,6,7],[5,7,8]]
    model = larExtrude((V,FV),2*[1,2,-3])
    VIEW(EXPLODE(1,1,1.2)(MKPOLS(model)))
    
    V0 = [[]]
    CV0 = [[0]]
    model = larExtrude((V0,CV0),6*[1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
    model = larExtrude(model,6*[1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
    model = larExtrude(model,6*[1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
    
    V0 = [[]]
    CV0 = [[0]]
    model = larExtrude((V0,CV0),2*[1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
    model = larExtrude(model,2*[1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))

#------------------------------------------------------------------
def simplexGrid(args):
    model = ([[]],[[0]])
    for k,steps in enumerate(args):
        model = larExtrude(model,steps*[1])
    V,cells = model
    verts = AA(list)(scipy.array(V) / AA(float)(args))
    return [verts, cells]

if __name__ == "__main__":
    
    grid_2d = simplexGrid([3,3])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(grid_2d)))
    
    grid_3d = simplexGrid([2,3,4])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(grid_3d)))

#------------------------------------------------------------------
def simplexFacets(simplices):
    """
        Estraction of non-oriented (d-1)-facets of
        d-dimensional "simplices".
        
        Return a list of d-tuples of integers
        """
    out = []
    d = len(simplices[0])+1
    for simplex in simplices:
        out += [simplex[0:k]+simplex[k+1:d] for k in range(d-1)]
    out = sorted(out)
    return [simplex for k,simplex in enumerate(out[:-1])
            if out[k] != out[k+1]] + [out[-1]]

if __name__ == "__main__":

    V,CV = simplexGrid([1,1,1])
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))
    SK2 = (V,simplexFacets(CV))
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK2)))
    SK1 = (V,simplexFacets(SK2[1]))
    VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK1)))
    
    print "\nk_0, k_1, k_2, k_3 =",len(V),len(SK1[1]),len(SK2[1]),len(CV), "\n"

#------------------------------------------------------------------

