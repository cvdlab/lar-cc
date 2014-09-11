
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from bool1 import *


n = 24
V1 = [[5*cos(angle*2*PI/n)+2.5, 5*sin(angle*2*PI/n)+2.5] for angle in range(n)]
FV1 = [range(n)]
EV1 = TRANS([range(n),range(1,n+1)]); EV1[-1] = [0,n-1]
VV1 = AA(LIST)(range(len(V1)))

V2 = [[4*cos(angle*2*PI/n), 4*sin(angle*2*PI/n)] for angle in range(n)]
FV2 = [range(n)]
EV2 = EV1
VV2 = AA(LIST)(range(len(V2)))


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

