""" Cylinder implemented as 9-point NURBS curve """
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/') 
from splines import *


knots = [0,0,0,1,1,2,2,3,3,4,4,4]
_p=math.sqrt(2)/2.0
controls = [[-1,0,1], [-_p,_p,_p], [0,1,1], [_p,_p,_p],[1,0,1], [_p,-_p,_p], 
         [0,-1,1], [-_p,-_p,_p], [-1,0,1]]
c1 = BEZIER(S1)([[-1,0,0,1],[-1,0,1,1]])
c2 = BEZIER(S1)([[-_p,_p,0,_p],[-_p,_p,_p,_p]])
c3 = BEZIER(S1)([[0,1,0,1],[0,1,1,1]])
c4 = BEZIER(S1)([[_p,_p,0,_p],[_p,_p,_p,_p]])
c5 = BEZIER(S1)([[1,0,0,1],[1,0,1,1]])
c6 = BEZIER(S1)([[_p,-_p,0,_p],[_p,-_p,_p,_p]])
c7 = BEZIER(S1)([[0,-1,0,1],[0,-1,1,1]])
c8 = BEZIER(S1)([[-_p,-_p,0,_p],[-_p,-_p,_p,_p]])
c9 = BEZIER(S1)([[-1,0,0,1],[-1,0,1,1]])
controls = [c1,c2,c3,c4,c5,c6,c7,c8,c9]
         
tnurbs = TNURBS(S2)(2)(knots)(controls)
dom = larModelProduct([larDomain([10]),larDom(knots)])
dom = larIntervals([10,36],'simplex')([1,4])
obj = larMap(tnurbs)(dom)
VIEW(STRUCT( MKPOLS(obj) ))
