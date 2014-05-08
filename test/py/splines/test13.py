""" Periodic B-spline curve """
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/') 
from splines import *


controls = [[0,1],[0,0],[1,0],[1,1],[0,1]]
knots = [0,0,0,1,2,3,3,3]           # non-periodic B-spline
tbspline = TBSPLINE(S1)(2)(knots)(controls)
obj = larMap(tbspline)(larDom(knots))  
VIEW(STRUCT( MKPOLS(obj) + [POLYLINE(controls)] ))

knots = [0,1,2,3,4,5,6,7]           # periodic B-spline
tbspline = TBSPLINE(S1)(2)(knots)(controls)  
obj = larMap(tbspline)(larDom(knots))
VIEW(STRUCT( MKPOLS(obj) + [POLYLINE(controls)] ))
""" Circle implemented as 9-point NURBS curve """
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/') 
from splines import *


knots = [0,0,0,1,1,2,2,3,3,4,4,4]
_p=math.sqrt(2)/2.0
controls = [[-1,0,1], [-_p,_p,_p], [0,1,1], [_p,_p,_p],[1,0,1], [_p,-_p,_p], 
         [0,-1,1], [-_p,-_p,_p], [-1,0,1]]
nurbs = NURBS(2)(knots)(controls)
obj = larMap(nurbs)(larDom(knots))
VIEW(STRUCT( MKPOLS(obj) + [POLYLINE(controls)] ))
