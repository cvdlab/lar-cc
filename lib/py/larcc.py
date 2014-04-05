# -*- coding: utf-8 -*-
""" Basic LARCC library """

"""
The MIT License
===============
    
Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
'Software'), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from pyplasm import *
import collections
import scipy
import numpy as np
from scipy import zeros,arange,mat,amin,amax
from scipy.sparse import vstack,hstack,csr_matrix,coo_matrix,lil_matrix,triu

from lar2psm import *

def triples2mat(triples,shape="csr"):
    n = len(triples)
    data = arange(n)
    ij = arange(2*n).reshape(2,n)
    for k,item in enumerate(triples):
        ij[0][k],ij[1][k],data[k] = item
    return scipy.sparse.coo_matrix((data, ij)).asformat(shape)

def brc2Coo(ListOfListOfInt):
    COOm = [[k,col,1] for k,row in enumerate(ListOfListOfInt)
            for col in row ]
    return COOm

def coo2Csr(COOm):
    CSRm = triples2mat(COOm,"csr")
    return CSRm

def csrCreate(BRCmatrix,shape=(0,0)):
    triples = brc2Coo(BRCmatrix)
    if shape == (0,0):
        CSRmatrix = coo2Csr(triples)
    else:
        CSRmatrix = scipy.sparse.csr_matrix(shape)
        for i,j,v in triples: CSRmatrix[i,j] = v
    return CSRmatrix

def csrGetNumberOfRows(CSRmatrix):
    Int = CSRmatrix.shape[0]
    return Int
    
def csrGetNumberOfColumns(CSRmatrix):
    Int = CSRmatrix.shape[1]
    return Int

def csr2DenseMatrix(CSRm):
    nrows = csrGetNumberOfRows(CSRm)
    ncolumns = csrGetNumberOfColumns(CSRm)
    ScipyMat = zeros((nrows,ncolumns),int)
    C = CSRm.tocoo()
    for triple in zip(C.row,C.col,C.data):
        ScipyMat[triple[0],triple[1]] = triple[2]
    return ScipyMat

def matrixProduct(CSRm1,CSRm2):
    CSRm = CSRm1 * CSRm2
    return CSRm

def csrTranspose(CSRm):
    CSRm = CSRm.T
    return CSRm

def csrBoundaryFilter(CSRm, facetLengths):
    maxs = [max(CSRm[k].data) for k in range(CSRm.shape[0])]
    inputShape = CSRm.shape
    coo = CSRm.tocoo()
    for k in range(len(coo.data)):
        if coo.data[k]==maxs[coo.row[k]]: coo.data[k] = 1
        else: coo.data[k] = 0
    mtx = coo_matrix((coo.data, (coo.row, coo.col)), shape=inputShape)
    out = mtx.tocsr()
    return out

def csrPredFilter(CSRm, pred):
   # can be done in parallel (by rows)
   coo = CSRm.tocoo()
   triples = [[row,col,val] for row,col,val 
            in zip(coo.row,coo.col,coo.data) if pred(val)]
   i, j, data = TRANS(triples)
   CSRm = scipy.sparse.coo_matrix((data,(i,j)),CSRm.shape).tocsr()
   return CSRm

def boundary(cells,facets):
    csrCV = csrCreate(cells)
    csrFV = csrCreate(facets)
    csrFC = matrixProduct(csrFV, csrTranspose(csrCV))
    facetLengths = [csrCell.getnnz() for csrCell in csrCV]
    return csrBoundaryFilter(csrFC,facetLengths)

def coboundary(cells,facets):
    Boundary = boundary(cells,facets)
    return csrTranspose(Boundary)

def zeroChain(cells):
   pass

def totalChain(cells):
   return csrCreate([[0] for cell in cells])

def boundaryCells(cells,facets):
    csrBoundaryMat = boundary(cells,facets)
    csrChain = totalChain(cells)
    csrBoundaryChain = matrixProduct(csrBoundaryMat, csrChain)
    for k,value in enumerate(csrBoundaryChain.data):
        if value % 2 == 0: csrBoundaryChain.data[k] = 0
    boundaryCells = [k for k,val in enumerate(csrBoundaryChain.data.tolist()) if val == 1]
    return boundaryCells

def signedBoundary (V,CV,FV):
   # compute the set of pairs of indices to [boundary face,incident coface]
   coo = boundary(CV,FV).tocoo()
   pairs = [[coo.row[k],coo.col[k]] for k,val in enumerate(coo.data) if val != 0]
   
   # compute the [face, coface] pair as vertex lists
   vertLists = [[FV[pair[0]], CV[pair[1]]]for pair in pairs]
   
   # compute two n-cells to compare for sign
   cellPairs = [ [list(set(coface).difference(face))+face,coface] 
               for face,coface in vertLists]
   
   # compute the local indices of missing boundary cofaces
   missingVertIndices = [ coface.index(list(set(coface).difference(face))[0]) 
                     for face,coface in vertLists]
   
   # compute the point matrices to compare for sign
   pointArrays = [ [[V[k]+[1.0] for k in facetCell], [V[k]+[1.0] for k in cofaceCell]] 
               for facetCell,cofaceCell in cellPairs]
   
   # signed incidence coefficients
   cofaceMats = TRANS(pointArrays)[1]
   cofaceSigns = AA(SIGN)(AA(np.linalg.det)(cofaceMats))
   faceSigns = AA(C(POWER)(-1))(missingVertIndices)
   signPairProd = AA(PROD)(TRANS([cofaceSigns,faceSigns]))
   
   # signed boundary matrix
   csrSignedBoundaryMat = csr_matrix( (signPairProd,TRANS(pairs)) )
   return csrSignedBoundaryMat

def signedBoundaryCells(verts,cells,facets):
   csrBoundaryMat = signedBoundary(verts,cells,facets)
   csrTotalChain = totalChain(cells)
   csrBoundaryChain = matrixProduct(csrBoundaryMat, csrTotalChain)
   coo = csrBoundaryChain.tocoo()
   boundaryCells = list(coo.row * coo.data)
   return AA(int)(boundaryCells)
def pivotSimplices(V,CV,d=3):
   simplices = []
   for cell in CV:
      vcell = np.array([V[v] for v in cell])
      m, simplex = len(cell), []
      # translate the cell: for each k, vcell[k] -= vcell[0], and simplex[0] := cell[0]
      for k in range(m-1,-1,-1): vcell[k] -= vcell[0]
      # simplex = [0], basis = [], tensor = Id(d+1)
      simplex += [cell[0]]
      basis = []
      tensor = np.array(IDNT(d))
      # look for most far cell vertex
      dists = [SUM([SQR(x) for x in v])**0.5 for v in vcell]
      maxDistIndex = max(enumerate(dists),key=lambda x: x[1])[0]
      vector = np.array([vcell[maxDistIndex]])
      # normalize vector
      den=(vector**2).sum(axis=-1) **0.5
      basis = [vector/den]
      simplex += [cell[maxDistIndex]]
      unUsedIndices = [h for h in cell if h not in simplex]
      
      # for k in {2,d+1}:
      for k in range(2,d+1):
         # update the orthonormal tensor
         e = basis[-1]
         tensor = tensor - np.dot(e.T, e)
         # compute the index h of a best vector
         # look for most far cell vertex
         dists = [SUM([SQR(x) for x in np.dot(tensor,v)])**0.5
         if h in unUsedIndices else 0.0
         for (h,v) in zip(cell,vcell)]
         # insert the best vector index h in output simplex
         maxDistIndex = max(enumerate(dists),key=lambda x: x[1])[0]
         vector = np.array([vcell[maxDistIndex]])
         # normalize vector
         den=(vector**2).sum(axis=-1) **0.5
         basis += [vector/den]
         simplex += [cell[maxDistIndex]]
         unUsedIndices = [h for h in cell if h not in simplex]
      simplices += [simplex]
   return simplices

def simplexOrientations(V,simplices):
   vcells = [[V[v]+[1.0] for v in simplex] for simplex in simplices]
   return [SIGN(np.linalg.det(vcell)) for vcell in vcells]

def larCellAdjacencies(CSRm):
    CSRm = matrixProduct(CSRm,csrTranspose(CSRm))
    return CSRm

def setup(model,dim):
    V, cells = model
    csr = csrCreate(cells)
    csrAdjSquareMat = larCellAdjacencies(csr)
    csrAdjSquareMat = csrPredFilter(csrAdjSquareMat, GE(dim)) # ? HOWTODO ?
    return V,cells,csr,csrAdjSquareMat

def larFacets(model,dim=3):
    """
        Estraction of (d-1)-cellFacets from "model" := (V,d-cells)
        Return (V, (d-1)-cellFacets)
      """
    V,cells,csr,csrAdjSquareMat = setup(model,dim)
    cellFacets = []
    # for each input cell i
    for i in range(len(cells)):
        adjCells = csrAdjSquareMat[i].tocoo()
        cell1 = csr[i].tocoo().col
        pairs = zip(adjCells.col,adjCells.data)
        for j,v in pairs:
            if (i<j):
                cell2 = csr[j].tocoo().col
                cell = list(set(cell1).intersection(cell2))
                cellFacets.append(sorted(cell))
    # sort and remove duplicates
    cellFacets = sorted(AA(list)(set(AA(tuple)(cellFacets))))
    return V,cellFacets


if __name__ == "__main__": 
   
   print "\n>>> brc2Coo"
   V = [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]]
   FV = [[0, 1, 3], [1, 2, 4], [1, 3, 4], [2, 4, 5]]
   EV = [[0,1],[0,3],[1,2],[1,3],[1,4],[2,4],[2,5],[3,4],[4,5]]
   cooFV = brc2Coo(FV)
   cooEV = brc2Coo(EV)
   assert cooFV == [[0,0,1],[0,1,1],[0,3,1],[1,1,1],[1,2,1],[1,4,1],[2,1,1],
   [2,3,1], [2,4,1],[3,2,1],[3,4,1],[3,5,1]]
   assert cooEV == [[0,0,1],[0,1,1],[1,0,1],[1,3,1],[2,1,1],[2,2,1],[3,1,1],
   [3,3,1],[4,1,1],[4,4,1],[5,2,1],[5,4,1],[6,2,1],[6,5,1],[7,3,1],[7,4,1],
   [8,4,1],[8,5,1]]
   
   csrFV = coo2Csr(cooFV)
   csrEV = coo2Csr(cooEV)
   print "\ncsr(FV) =\n", repr(csrFV)
   print "\ncsr(EV) =\n", repr(csrEV)
   
   print "\n>>> brc2Csr"
   V = [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]]
   FV = [[0, 1, 3], [1, 2, 4], [1, 3, 4], [2, 4, 5]]
   EV = [[0,1],[0,3],[1,2],[1,3],[1,4],[2,4],[2,5],[3,4],[4,5]]
   csrFV = csrCreate(FV)
   csrEV = csrCreate(EV)
   print "\ncsrCreate(FV) =\n", csrFV
   VIEW(STRUCT(MKPOLS((V,FV))))
   VIEW(STRUCT(MKPOLS((V,EV))))
   
   print "\n>>> csrGetNumberOfRows"
   print "\ncsrGetNumberOfRows(csrFV) =", csrGetNumberOfRows(csrFV)
   print "\ncsrGetNumberOfRows(csrEV) =", csrGetNumberOfRows(csrEV)
   print "\n>>> csrGetNumberOfColumns"
   print "\ncsrGetNumberOfColumns(csrFV) =", csrGetNumberOfColumns(csrFV)
   print "\ncsrGetNumberOfColumns(csrEV) =", csrGetNumberOfColumns(csrEV)
   
   print "\n>>> csr2DenseMatrix"
   print "\nFV =\n", csr2DenseMatrix(csrFV)
   print "\nEV =\n", csr2DenseMatrix(csrEV)
   
   print "\n>>> csrBoundaryFilter"
   csrEF = matrixProduct(csrFV, csrTranspose(csrEV)).T
   facetLengths = [csrCell.getnnz() for csrCell in csrEV]
   CSRm = csrBoundaryFilter(csrEF, facetLengths).T
   print "\ncsrMaxFilter(csrFE) =\n", csr2DenseMatrix(CSRm)
   
   print "\n>>> csrPredFilter"
   CSRm = csrPredFilter(matrixProduct(csrFV, csrTranspose(csrEV)).T, GE(2)).T
   print "\nccsrPredFilter(csrFE) =\n", csr2DenseMatrix(CSRm)
   
   V = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], 
   [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
   
   CV =[[0, 1, 2, 4], [1, 2, 4, 5], [2, 4, 5, 6], [1, 2, 3, 5], [2, 3, 5, 6], 
   [3, 5, 6, 7]]
   
   FV =[[0, 1, 2], [0, 1, 4], [0, 2, 4], [1, 2, 3], [1, 2, 4], [1, 2, 5], 
   [1, 3, 5], [1, 4, 5], [2, 3, 5], [2, 3, 6], [2, 4, 5], [2, 4, 6], [2, 5, 6], 
   [3, 5, 6], [3, 5, 7], [3, 6, 7], [4, 5, 6], [5, 6, 7]]
   
   EV =[[0, 1], [0, 2], [0, 4], [1, 2], [1, 3], [1, 4], [1, 5], [2, 3], [2, 4], 
   [2, 5], [2, 6], [3, 5], [3, 6], [3, 7], [4, 5], [4, 6], [5, 6], [5, 7], 
   [6, 7]]
   
   print "\ncoboundary_2 =\n", csr2DenseMatrix(coboundary(CV,FV))
   print "\ncoboundary_1 =\n", csr2DenseMatrix(coboundary(FV,EV))
   print "\ncoboundary_0 =\n", csr2DenseMatrix(coboundary(EV,AA(LIST)(range(len(V)))))
   
   boundaryCells_2 = boundaryCells(CV,FV)
   boundaryCells_1 = boundaryCells([FV[k] for k in boundaryCells_2],EV)
   
   print "\nboundaryCells_2 =\n", boundaryCells_2
   print "\nboundaryCells_1 =\n", boundaryCells_1
   
   boundary = (V,[FV[k] for k in boundaryCells_2])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(boundary)))
   
   print "\n>>> larCellAdjacencies"
   adj_2_cells = larCellAdjacencies(csrFV)
   print "\nadj_2_cells =\n", csr2DenseMatrix(adj_2_cells)
   adj_1_cells = larCellAdjacencies(csrEV)
   print "\nadj_1_cells =\n", csr2DenseMatrix(adj_1_cells)
   
   V = [[0.,0.],[3.,0.],[0.,3.],[3.,3.],[1.,2.],[2.,2.],[1.,1.],[2.,1.]]
   FV = [[0,1,6,7],[0,2,4,6],[4,5,6,7],[1,3,5,7],[2,3,4,5],[0,1,2,3]]
   
   _,EV = larFacets((V,FV),dim=2)
   print "\nEV =",EV
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,EV))))
   
   FV = [[0,1,3],[1,2,4],[2,4,5],[3,4,6],[4,6,7],[5,7,8], # full
      [1,3,4],[4,5,7], # empty
      [0,1,2],[6,7,8],[0,3,6],[2,5,8]] # exterior
         
   _,EV = larFacets((V,FV),dim=2)
   print "\nEV =",EV
   
   
