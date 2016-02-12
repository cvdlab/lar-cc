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
def csrBoundaryFilter1(csrBoundaryBoundaryMat,cells,facets,faces,lenV,FE, CSRm, facetLengths):
    maxs = [max(CSRm[k].data) for k in range(CSRm.shape[0])]
    inputShape = CSRm.shape
    coo = CSRm.tocoo()
    for k in range(len(coo.data)):
        if coo.data[k]==maxs[coo.row[k]]: coo.data[k] = 1
        else: coo.data[k] = 0
    mtx = coo_matrix((coo.data, (coo.row, coo.col)), shape=inputShape)
    out = mtx.tocsr()
    
    unreliable = [k for k in range(out.shape[0]) if sum(out[k,:].todense()[0]) > 2]
    if unreliable != []:
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
    csrCV = csrCreate(CV,lenV)
    csrFV = csrCreate(FV,lenV)
    csrFC = csrFV * csrCV.T
    facetLengths = [csrCell.getnnz() for csrCell in csrCV]
    VV = AA(LIST)(range(lenV))
    csrBBMat = boundary(FV,EV)*boundary(CV,FV)
    FE = larcc.crossRelation(lenV,FV,EV,True)
    return csrBoundaryFilter1(csrBBMat,CV,FV,EV,lenV,FE,csrFC,facetLengths)

