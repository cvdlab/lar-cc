""" Skeletons and oriented boundary of a simplicial complex """
from larlib import *


""" Oriented boundary matrix visualization """
np.set_printoptions(threshold='nan')
csrSignedBoundaryMat = signedSimplicialBoundary (V,FV,EV)
Z = csr2DenseMatrix(csrSignedBoundaryMat)
print "\ncsrSignedBoundaryMat =\n", Z
import matplotlib.pyplot
from pylab import *
matshow(Z)
show()

"""  Computation of oriented boundary cells """
boundaryCells_1 = signedBoundaryCells(V,FV,EV)
print "\nboundaryCells_1 =\n", boundaryCells_1
boundaryEV = [EV[-k] if k<0 else swap(EV[k]) for k in boundaryCells_1]
bndry = (V,boundaryEV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(bndry)))

