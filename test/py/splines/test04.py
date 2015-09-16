""" Graph of Bernstein-Bezier basis """
from larlib import *

def larBezierBasisGraph(degree):
   basis = larBernsteinBasis(S1)(degree)
   dom = larDomain([32])
   graphs = CONS(AA(larMap)(DISTL([S1, basis])))(dom)
   return graphs

graphs = larBezierBasisGraph(4)
VIEW(STRUCT( CAT(AA(MKPOLS)( graphs )) ))
