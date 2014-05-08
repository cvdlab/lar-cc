""" Effect of knot multiplicity on B-spline curve """
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/') 
from splines import *


points = [[0,0],[-1,2],[1,4],[2,3],[1,1],[1,2],[2.5,1]]
b1 = BSPLINE(2)([0,0,0,1,2,3,4,5,5,5])(points)
VIEW(STRUCT(MKPOLS( larMap(b1)(larDom([0,5])) ) + [POLYLINE(points)]))
b2 = BSPLINE(2)([0,0,0,1,1,2,3,4,4,4])(points)
VIEW(STRUCT(MKPOLS( larMap(b2)(larDom([0,5])) ) + [POLYLINE(points)]))
b3 = BSPLINE(2)([0,0,0,1,1,1,2,3,3,3])(points)
VIEW(STRUCT(MKPOLS( larMap(b3)(larDom([0,5])) ) + [POLYLINE(points)]))
b4 = BSPLINE(2)([0,0,0,1,1,1,1,2,2,2])(points)
VIEW(STRUCT(MKPOLS( larMap(b4)(larDom([0,5])) ) + [POLYLINE(points)]))
