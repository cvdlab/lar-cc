
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from bool import *
""" Definition of Boolean arguments """
V1 = [[3,0],[11,0], [13,10], [10,11], [8,11], [6,11], [4,11], [1,10], [4,3], [6,4], 
   [8,4], [10,3]]
FV1 = [[0,1,8,9,10,11],[1,2,11], [3,10,11], [4,5,9,10], [6,8,9], [0,7,8], [2,3,11],
   [3,4,10], [5,6,9], [6,7,8], range(8)]
V2 = [[0,3],[14,2], [14,5], [14,7], [14,11], [0,8], [3,7], [3,5]]
FV2 =[[0,5,6,7], [0,1,7], [4,5,6], [2,3,6,7], [1,2,7], [3,4,6], range(6)]

""" Bulk of Boolean task computation """
""" Computation of edges an input visualisation """
model1 = V1,FV1
model2 = V2,FV2
submodel = SKEL_1(STRUCT(MKPOLS(model1)+MKPOLS(model2)))
VV1 = AA(LIST)(range(len(V1)))
_,EV1 = larFacets((V1,FV1),dim=2,emptyCellNumber=1)
VV2 = AA(LIST)(range(len(V2)))
_,EV2 = larFacets((V2,FV2),dim=2,emptyCellNumber=1)
VIEW(larModelNumbering(V1,[VV1,EV1,FV1],submodel,4))
VIEW(larModelNumbering(V2,[VV2,EV2,FV2],submodel,4))

V,[VV,EEV1,EEV2,CV1,CV2],n12 = covering(model1,model2)
CCV = CV1+CV2
EEV = EEV1+EEV2
VIEW(larModelNumbering(V,[VV,EEV,CCV],submodel,4))

CV, BV1, BV2, BF1, BF2, BV, BF, nE1 = partition(V, CV1,CV2, EEV1,EEV2)
boundaries = COLOR(YELLOW)(SKEL_1(STRUCT(MKPOLS((V,[EEV[e] for e in BF])))))
submodel = STRUCT([ SKEL_1(STRUCT(MKPOLS((V,CV)))), boundaries ])
VIEW(larModelNumbering(V,[VV,EEV,CV],submodel,4))
""" Inversion of incidences """
VC = invertRelation(V,CV)
VC1 = invertRelation(V,CV1)
VC2 = invertRelation(V,CV2)
VEE1 = invertRelation(V,EEV1)
VEE2 = [[e+nE1  for e in vE] for vE in invertRelation(V,EEV2)]
submodel = SKEL_1(STRUCT(MKPOLS((V,CV1+CV2))))
VE = [VEE1[v]+VEE2[v] for v in range(len(V))]


n0,n1 = 0, max(AA(max)(CV1))        # vertices in CV1 (extremes included)
m0,m1 = n1+1-n12, max(AA(max)(CV2))    # vertices in CV2 (extremes included)
VE = [VEE1[v]+VEE2[v] for v in range(len(V))]
cells12 = mixedCells(CV,n0,n1,m0,m1)
pivots = mixedCellsOnBoundaries(cells12,BV1,BV2)
tasks = splittingTasks(V,pivots,BV,BF,VC,CV,EEV,VE)
print "\ntasks (face,cell) =",tasks
   
   
dict_fc,dict_cf = initTasks(tasks)
print "\ntasks (dict_fc) =",dict_fc
print "\ntasks (dict_cf) =",dict_cf

vertdict,cellPairs,nverts = splitCellsCreateVertices(dict_fc,dict_cf,V,EEV,CV,VC)

V = list(None for k in range(nverts))  
for item in vertdict.items(): V[item[1][0]] = eval(item[0])
VV = AA(LIST)(range(len(V)))
cells1,cells2 = TRANS(cellPairs)
out = [COLOR(WHITE)(MKPOL([V,[[v+1 for v in cell] for cell in cells1],None])), 
      COLOR(MAGENTA)(MKPOL([V,[[v+1 for v in cell] for cell in cells2],None]))]

boundaries = COLOR(YELLOW)(SKEL_1(STRUCT(MKPOLS((V,[EEV[e] for e in BF])))))
submodel = STRUCT([ SKEL_1(STRUCT(MKPOLS((V,CV)))), boundaries ])
VIEW(STRUCT([ STRUCT(out), larModelNumbering(V,[VV,[],CV],submodel,4) ]))

