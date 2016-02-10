""" boundary operators """
from larlib import *
""" convex-cells boundary operator """
from larcc import *
def boundary(cells,facets):
    lenV = max(CAT(cells))+1
    csrCV = csrCreate(cells,lenV)
    csrFV = csrCreate(facets,lenV)
    csrFC = matrixProduct(csrFV, csrTranspose(csrCV))
    facetLengths = [csrCell.getnnz() for csrCell in csrCV]
    return csrBoundaryFilter(csrFC,facetLengths)

""" path-connected-cells boundary operator """
import larcc
from larcc import *
def boundary1(CV,FV,EV):
    lenV = max(CAT(CV))+1
    csrCV = csrCreate(CV,lenV)
    csrFV = csrCreate(FV,lenV)
    csrFC = matrixProduct(csrFV, csrTranspose(csrCV))
    facetLengths = [csrCell.getnnz() for csrCell in csrCV]
    VV = AA(LIST)(range(lenV))
    csrBBMat = boundary(FV,EV)*boundary(CV,FV)
    FE = larcc.crossRelation(lenV,FV,EV,True)
    return csrBoundaryFilter1(csrBBMat,CV,FV,EV,lenV,FE,csrFC,facetLengths)

