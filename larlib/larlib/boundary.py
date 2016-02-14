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

def csrBoundaryFilter1(csrBoundaryBoundaryMat,cells,facets,faces,lenV,FE):
    out = boundary(cells,facets)
    
    def csrRowSum(h): return sum(out.data[out.indptr[h]:out.indptr[h+1]])
    
    unreliable = [h for h in range(out.shape[0]) if csrRowSum(h) > 2]
    if unreliable != []:
        print "\n>>>>> unreliable =",unreliable
        for row in unreliable:
            for j in range(len(cells)):
                if out[row,j] == 1:
                    csrCFE = csrBoundaryBoundaryMat[:,j]
                    cooCFE = csrCFE.tocoo()
                    flawedCells = [cooCFE.row[k] for k,datum in enumerate(cooCFE.data)
                        if datum>2]
                    if all([facet in flawedCells  for facet in FE[row]]):
                        out[row,j]=0
    return out

def boundary1(CV,FV,EV):
    lenV = max(CAT(CV))+1
    csrBBMat = boundary(FV,EV) * boundary(CV,FV)
    print "\ncsrBBMat =",csrBBMat,"\n"
    FE = larcc.crossRelation(lenV,FV,EV)
    return csrBoundaryFilter1(csrBBMat,CV,FV,EV,lenV,FE)

