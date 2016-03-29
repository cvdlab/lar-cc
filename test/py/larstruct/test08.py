""" LAR model input and handling """
from larlib import *

""" Input of LAR architectural plan """
from larlib import *

V = [[3,-3],
[9,-3],[0,0],[3,0],[9,0],[15,0],
[3,3],[6,3],[9,3],[15,3],[21,3], 
[0,9],[6,9],[15,9],[18,9],[0,13],
[6,13],[9,13],[15,13],[18,10],[21,10], 
[18,13],[6,16],[9,16],[9,17],[15,17],
[18,17],[-3,24],[6,24],[15,24],[-3,13]]
FV = [
[22,23,24,25,29,28], [15,16,22,28,27,30], [18,21,26,25], 
[13,14,19,21,18], [16,17,23,22], [11,12,16,15],
[9,10,20,19,14,13], [2,3,6,7,12,11], [0,1,4,8,7,6,3],
[4,5,9,13,18,17,16,12,7,8],[17,18,25,24,23]]

polylines = lar2polylines((V,FV))
lines = CAT([zip(polyline[:-1],polyline[1:]) for polyline in polylines])   
verts = dict(zip(AA(vcode(4))(V),range(len(V))))
edges = [tuple(sorted([verts[vcode(4)(v1)], verts[vcode(4)(v2)]])) for v1,v2 in lines]
EV = list(set(edges))

dwelling = larApply(t(3,0))(Model((V,FV)))
print "\n dwelling =",dwelling
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((dwelling.verts,dwelling.cells))))
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((dwelling.verts,EV))))
plan = Struct([dwelling,s(-1,1),dwelling])
VIEW(EXPLODE(1.2,1.2,1)(CAT(AA(MKPOLS)(evalStruct(plan)))))
