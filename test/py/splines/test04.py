""" Graph of Bernstein-Bezier basis """
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/') 
from splines import *


def larBezierBasisGraph(degree):
   basis = larBernsteinBasis(S1)(degree)
   dom = larDomain([32])
   graphs = CONS(AA(larMap)(DISTL([S1, basis])))(dom)
   return graphs

graphs = larBezierBasisGraph(4)
VIEW(STRUCT( CAT(AA(MKPOLS)( graphs )) ))
