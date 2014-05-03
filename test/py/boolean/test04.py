""" test program for the boolean module """
from pyplasm import *
from scipy import *
import os,sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from lar2psm import *
from simplexn import *
from larcc import *
from largrid import *
from myfont import *
from mapper import *

from boolean import *
blue = larSimplexGrid([30,60])
V2,CV2 = larSimplexGrid([70,40])
V2 = larTranslate( [.5,.5])(V2)
red = V2,CV2
VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(blue) ))
VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(red) ))
V, CV1, CV2, n12 = vertexSieve(red,blue)
V, n1,n2,n12 = boolOps(red,blue)
CV = Delaunay(array(V)).vertices
