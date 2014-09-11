
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from bool1 import *

V1 = [[0,0,0],[10,0,0],[10,10,0],[0,10,0],[0,0,10],[10,0,10],[10,10,10],[0,10,10]]
V1,[VV1,EV1,FV1,CV1] = larCuboids((1,1,1),True)
V1 = [SCALARVECTPROD([5,v]) for v in V1]

V2 = [SUM([v,[2.5,2.5,2.5]]) for v in V1]
[VV2,EV2,FV2,CV2] = [VV1,EV1,FV1,CV1]


if DEBUG: VIEW(STRUCT(MKPOLS((V1,EV1)) + MKPOLS((V2,EV2))))

model1,model2 = (V1,FV1),(V2,FV2)
V, CV1,CV2, n1,n12,n2 = mergeVertices(model1,model2)  #<<<<<<<<<<<<<<<<

print "\nV =", V
print "\nCV1,CV2 =", CV1,CV2
print "\nn1,n12,n2 =", n1,n12,n2

submodel = SKEL_1(STRUCT(MKPOLS((V,CV1+CV2)))) 
VV = AA(LIST)(range(len(V)))
if DEBUG: VIEW(STRUCT([ submodel,larModelNumbering(V,[VV,_,CV1+CV2],submodel,3)]))

V,CV,vertDict,n1,n12,n2 = makeCDC(model1, model2)     #<<<<<<<<<<<<<<<<

print "\nCV =", CV
print "\nvertDict =", vertDict

submodel = SKEL_1(STRUCT(MKPOLS((V,CV))))
VIEW(STRUCT([ submodel,larModelNumbering(V,[VV,_,CV],submodel,4)]))

