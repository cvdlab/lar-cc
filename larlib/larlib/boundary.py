""" boundary operators """
from larlib import *
""" convex-cells boundary operator --- best implementation """
def boundary(cells,facets):
    lenV = max(CAT(cells))+1
    csrCV = csrCreate(cells,lenV)
    csrFV = csrCreate(facets,lenV)
    csrFC = csrFV * csrCV.T
    facetLengths = [csrFacet.getnnz() for csrFacet in csrFV]
    m,n = csrFC.shape
    facetCoboundary = [[csrFC.indices[csrFC.indptr[h]+k] 
        for k,v in enumerate(csrFC.data[csrFC.indptr[h]:csrFC.indptr[h+1]]) 
            if v==facetLengths[h]] for h in range(m)]
    indptr = [0]+list(cumsum(AA(len)(facetCoboundary)))
    indices = CAT(facetCoboundary)
    data = [1]*len(indices)
    return csr_matrix((data,indices,indptr),shape=(m,n),dtype='b')

""" path-connected-cells boundary operator """
import larlib
import larcc
from larcc import *

def csrBoundaryFilter1(unreliable,out,csrBBMat,cells,FE):
    for row in unreliable:
        for j in range(len(cells)):
            if out[row,j] == 1:
                cooCE = csrBBMat.T[j].tocoo()
                flawedCells = [cooCE.col[k] for k,datum in enumerate(cooCE.data)
                    if datum>2]
                if all([facet in flawedCells  for facet in FE[row]]):
                    out[row,j]=0
    return out

def boundary1(CV,FV,EV):
    def csrRowSum(h): 
        return sum(out.data[out.indptr[h]:out.indptr[h+1]])    
    lenV = max(CAT(CV))+1
    csrBBMat = boundary(FV,EV) * boundary(CV,FV)
    print "\ncsrBBMat =",csrBBMat,"\n"
    FE = larcc.crossRelation(lenV,FV,EV)
    out = boundary(CV,FV)
    unreliable = [h for h in range(len(FV)) if csrRowSum(h) > 2]
    if unreliable != []:
        out = csrBoundaryFilter1(unreliable,out,csrBBMat,CV,FE)
    return out

