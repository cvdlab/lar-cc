""" Periodic B-spline curve """
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/') 
from splines import *


controls = [[0,1],[0,0],[1,0],[1,1],[0,1]]
knots = [0,0,0,1,2,3,3,3]           # non-periodic B-spline
bspline = BSPLINE(2)(knots)(controls)
obj = larMap(bspline)(larDom(knots))  
VIEW(STRUCT( MKPOLS(obj) + [POLYLINE(controls)] ))

knots = [0,1,2,3,4,5,6,7]           # periodic B-spline
bspline = BSPLINE(2)(knots)(controls)  
obj = larMap(bspline)(larDom(knots))
VIEW(STRUCT( MKPOLS(obj) + [POLYLINE(controls)] ))
