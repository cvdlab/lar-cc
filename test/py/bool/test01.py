import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from bool import *
""" Definition of Boolean arguments """
V1 = [[3,0],[11,0], [13,10], [10,11], [8,11], [6,11], [4,11], [1,10], [4,3], [6,4], 
   [8,4], [10,3]]
FV1 = [[0,1,8,9,10,11],[1,2,11], [3,10,11], [4,5,9,10], [6,8,9], [0,7,8], [2,3,11],
   [3,4,10], [5,6,9], [6,7,8]]
EV1 = [[0,1],[0,7],[0,8],[1,2],[1,11],[2,3],[2,11],[3,4],[3,10],[3,11],[4,5],[4,10],[5,6],[5,9],[6,7],[6,8],[6,9],[7,8],[8,9],[9,10],[10,11]]
VV1 = AA(LIST)(range(len(V1)))

V2 = [[0,3],[14,2], [14,5], [14,7], [14,11], [0,8], [3,7], [3,5]]
FV2 =[[0,5,6,7], [0,1,7], [4,5,6], [2,3,6,7], [1,2,7], [3,4,6]]
EV2 = [[0,1],[0,5],[0,7],[1,2],[1,7],[2,3],[2,7],[3,4],[3,6],[4,5],[4,6],[5,6],[6,7]]
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

BC1 = signedCellularBoundaryCells(V1,basis1)
print "\nsignedBoundaryCells1 =",BC1
BC2 = signedCellularBoundaryCells(V2,basis2)
print "\nsignedBoundaryCells2 =",BC2
BC = [[ vertdict[vcode(V1[v])][0] for v in cell] for cell in BC1] + [ [ vertdict[vcode(V2[v])][0] for v in cell] for cell in BC2]
BC = sorted(BC)

BV1 = list(set(CAT(BC1)))
BV1 = [vertdict[vcode(V1[v])][0] for v in BV1]
BV2 = list(set(CAT(BC2)))
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

cellPairs,twoCellIndices,cuttingFaces = splitCellsCreateVertices(vertdict,dict_fc,dict_cf,V,BC,CV,VC)
showSplitting(V,cellPairs,BC,CV)

splitCellsBits(cuttingFaces,cellPairs,twoCellIndices,CV1,CV2,n12,BC)
