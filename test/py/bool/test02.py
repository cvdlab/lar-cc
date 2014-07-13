from pyplasm import *
from scipy import *
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from lar2psm import *
from simplexn import *
from larcc import *
from largrid import *
from myfont import *
from mapper import *

from bool import *
V1 = [[1,1],[3,3],[3,1],[2,3],[2,1],[1,3]]
V2 = [[1,1],[1,3],[2,3],[2,2],[3,2],[0,1],[0,0],[2,0],[3,0]]
CV1 = [[0,3,4,5],[1,2,3,4]]
CV2 = [[3,4,7,8],[0,1,2,3,5,6,7]]
model1 = V1,CV1; model2 = V2,CV2
VIEW(STRUCT([ 
   COLOR(CYAN)(SKEL_1(STRUCT(MKPOLS(model1)))), 
   COLOR(RED)(SKEL_1(STRUCT(MKPOLS(model2)))) ]))
# V, n1,n2,n12,BV1,BV2 = boolOps(model1,model2)
# VIEW(SKEL_1(STRUCT(MKPOLS((V, CV_un[:n1]+CV_int )))))
# VIEW(SKEL_1(STRUCT(MKPOLS((V, CV_un[n1-n12:]+CV_int )))))
""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')
from bool import *

V1 = [[0,0],[10,0],[10,10],[0,10]]
FV1 = [range(4)]
EV1 = [[0,1],[1,2],[2,3],[3,0]]
VV1 = AA(LIST)(range(len(V1)))

V2 = [[2.5,2.5],[12.5,2.5],[12.5,12.5],[2.5,12.5]]
FV2 = [range(4)]
EV2 = [[0,1],[1,2],[2,3],[3,0]]
VV2 = AA(LIST)(range(len(V2)))
""" Bulk of Boolean task computation """
""" Computation of edges an input visualisation """
model1 = V1,FV1
model2 = V2,FV2
basis1 = [VV1,EV1,FV1]
basis2 = [VV2,EV2,FV2]
submodel12 = STRUCT(MKPOLS((V1,EV1))+MKPOLS((V2,basis2[1])))
VIEW(larModelNumbering(V1,basis1,submodel12,4))
VIEW(larModelNumbering(V2,basis2,submodel12,4))


V,[VV,_,_,CV1,CV2],n12 = covering(model1,model2,2,0)
CV = sorted(AA(sorted)(Delaunay(array(V)).vertices))
vertdict = defaultdict(list)
for k,v in enumerate(V): vertdict[vcode(v)] += [k]

BC1 = boundaryCells(basis1[-1],basis1[-2])
BC2 = boundaryCells(basis2[-1],basis2[-2])
BC = [[ vertdict[vcode(V1[v])][0] for v in basis1[-2][cell]] for cell in BC1] + [ [ vertdict[vcode(V2[v])][0] for v in basis2[-2][cell]] for cell in BC2]
BC = sorted(AA(sorted)(BC))

BV1 = list(set(CAT([basis1[-2][bc] for bc in BC1])))
BV1 = [vertdict[vcode(V1[v])][0] for v in BV1]
BV2 = list(set(CAT([basis2[-2][bc] for bc in BC2])))
BV1 = [vertdict[vcode(V2[v])][0] for v in BV2]
BV = list(set(CAT([v for v in BC])))
VV = AA(LIST)(range(len(V)))
submodel = STRUCT([SKEL_1(STRUCT(MKPOLS((V,CV)))), COLOR(RED)(STRUCT(MKPOLS((V,BC))))])
VIEW(larModelNumbering(V,[VV,BC,CV],submodel,4))

cells12 = mixedCells(CV,CV1,CV2,n12)
pivots = mixedCellsOnBoundaries(cells12,BV)
VBC = invertRelation(V,BC)
VC = invertRelation(V,CV)
tasks = splittingTasks(V,pivots,BV,BC,VBC,CV,VC)
dict_fc,dict_cf = initTasks(tasks)

cellPairs = splitCellsCreateVertices(vertdict,dict_fc,dict_cf,V,BC,CV,VC)
showSplitting(V,cellPairs,BC,CV)

