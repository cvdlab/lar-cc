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
import larcc
from larcc import *
def boundary1(CV,FV,EV):
    lenV = max(CAT(CV))+1
    csrCV = csrCreate(CV,lenV)
    csrFV = csrCreate(FV,lenV)
    csrFC = csrFV * csrCV.T
    facetLengths = [csrCell.getnnz() for csrCell in csrCV]
    VV = AA(LIST)(range(lenV))
    csrBBMat = boundary(FV,EV)*boundary(CV,FV)
    FE = larcc.crossRelation(lenV,FV,EV,True)
    return csrBoundaryFilter1(csrBBMat,CV,FV,EV,lenV,FE,csrFC,facetLengths)

