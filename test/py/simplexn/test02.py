
""" import modules from larcc/lib """
import sys
from scipy import reshape
sys.path.insert(0, 'lib/py/')
from pyplasm import *
from lar2psm import *
from simplexn import *

V,CV = larSimplexGrid1([2,2,2])
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))


FV = larSimplexFacets(CV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))
EV = larSimplexFacets(FV)
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,EV))))

