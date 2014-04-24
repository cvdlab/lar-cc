
""" import modules from larcc/lib """
import sys
from scipy import reshape
sys.path.insert(0, 'lib/py/')
from pyplasm import *
from lar2psm import *
from simplexn import *

V,FV = larSimplexGrid1([3,3])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))

EV = larSimplexFacets(FV)
ex,ey,ez = 1.5,1.5,1.5
VIEW(EXPLODE(ex,ey,ez)(MKPOLS((V,EV))))

