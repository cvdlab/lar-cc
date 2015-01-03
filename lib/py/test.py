    from pyplasm import *
from scipy import *
import os,sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from boolean import *
from matrix import *

subspace = T([1,2,3])([-50,-50,0])(CUBOID([100,100,100]))
VIEW(subspace)
#--------------------------------------------------------------------------

V1 = [[3,0],[11,0], [13,10], [10,11], [8,11], [6,11], [4,11], [1,10], [4,3], [6,4], 
	[8,4], [10,3]]
FV1 = [[0,1,8,9,10,11],[1,2,11], [3,10,11], [4,5,9,10], [6,8,9], [0,7,8], [2,3,11],
	[3,4,10], [5,6,9], [6,7,8], range(8)]
V2 = [[0,3],[14,2], [14,5], [14,7], [14,11], [0,8], [3,7], [3,5]]
FV2 =[[0,5,6,7], [0,1,7], [4,5,6], [2,3,6,7], [1,2,7], [3,4,6], range(6)]
model1 = V1,FV1
model2 = V2,FV2
submodel = SKEL_1(STRUCT(MKPOLS(model1)+MKPOLS(model2)))

VV1 = AA(LIST)(range(len(V1)))
_,EV1 = larFacets((V1,FV1),dim=2,emptyCellNumber=1)
VV2 = AA(LIST)(range(len(V2)))
_,EV2 = larFacets((V2,FV2),dim=2,emptyCellNumber=1)
"""
VIEW(SKEL_1(STRUCT(MKPOLS(model1))))
VIEW(SKEL_1(STRUCT(MKPOLS(model2))))
VIEW(submodel)
VIEW(larModelNumbering(V1,[VV1,EV1,FV1],submodel,4))
VIEW(larModelNumbering(V2,[VV2,EV2,FV2],submodel,4))
"""
V, CV1, CV2, n12 = vertexSieve(model1,model2)
_,EEV1 = larFacets((V,CV1),dim=2,emptyCellNumber=1)
_,EEV2 = larFacets((V,CV2),dim=2,emptyCellNumber=1)
CV1 = CV1[:-1]
CV2 = CV2[:-1]
VV = AA(LIST)(range(len(V)))
EEV = EEV1+EEV2
VIEW(larModelNumbering(V,[VV,EEV1+EEV2,CV1+CV2],submodel,4))

CV = sorted(AA(sorted)(Delaunay(array(V)).vertices))
CVdict = dict(zip(AA(tuple)(CV),range(len(CV))))
BV1, BV2,BF1,BF2 = boundaryVertices( V, CV1,CV2, 'cuboid', EEV1,EEV2 )
BV = BV1+BV2
"""
VIEW(STRUCT([
	EXPLODE(1.2,1.2,1)(MKPOLS((V,CV))),
	T(3)(.1)(COLOR(BLUE)(SKEL_1(EXPLODE(1.2,1.2,1)(MKPOLS((V,CV1)))))),
	T(3)(.1)(COLOR(MAGENTA)(SKEL_1(EXPLODE(1.2,1.2,1)(MKPOLS((V,CV2)))))),
	COLOR(WHITE)(EXPLODE(1.2,1.2,1)(AA(MK)([V[v] for v in BV1]))),
	COLOR(RED)(EXPLODE(1.2,1.2,1)(AA(MK)([V[v] for v in BV2])))
]))
"""
submodel = SKEL_1(STRUCT(MKPOLS((V,CV))))
VIEW(larModelNumbering(V,[VV,[],CV],submodel,4))

# transposition (relational inversion) of CV,CV1,CV2
def invertRelation(V,CV):
	VC = [[] for k in range(len(V))]
	for k,cell in enumerate(CV):
		for v in cell:
			VC[v] += [k]
	return VC
VC = invertRelation(V,CV)
VC1 = invertRelation(V,CV1)
VC2 = invertRelation(V,CV2)
submodel = SKEL_1(STRUCT(MKPOLS((V,CV1+CV2))))
VIEW(larModelNumbering(V,[VV,[],CV1],submodel,4))
VIEW(larModelNumbering(V,[VV,[],CV2],submodel,4))
VEE1 = invertRelation(V,EEV1)
nE1 = len(EV1)
VEE2 = [[e+nE1  for e in vE] for vE in invertRelation(V,EEV2)]

n0,n1 = 0, max(AA(max)(CV1))			# vertices in CV1 (extremes included)
m0,m1 = n1+1-n12, max(AA(max)(CV2))		# vertices in CV2 (extremes included)
VE = [VEE1[v]+VEE2[v] for v in range(len(V))]
VIEW(larModelNumbering(V,[VV,EEV1+EEV2,CV2],submodel,4))
VIEW(larModelNumbering(V,[VV,EEV1+EEV2,CV1+CV2],submodel,4))

# Constructing the queue of vertex-based jobs
"""
for v in BV:
	print VC[v],VC1[v],VC2[v]
	VIEW(EXPLODE(1.2,1.2,1)(
		AA(SKEL_1)(MKPOLS((V,CV))) + MKPOLS((V,[CV[cell] for cell in VC[v]])) +
		AA(COMP([COLOR(BLUE),SKEL_1]))(MKPOLS((V,[CV1[cell] for cell in VC1[v]]))) + 
		AA(COMP([COLOR(RED),SKEL_1]))(MKPOLS((V,[CV2[cell] for cell in VC2[v]])))
	))
"""
# Incidence on boundary cocells
cells = CV1; facets = EEV1
boundaryCells1, boundaryCocells1 = boundaryCellsCocells(cells,facets)
print "\nboundaryCocells1 =",boundaryCocells1
cells = CV2; facets = EEV2
boundaryCells2, boundaryCocells2 = boundaryCellsCocells(cells,facets)
print "\nboundaryCocells2 =",boundaryCocells2

# 1. Look for cells in Delaunay, with vertices in both operands
cells12 = [list(cell) for cell in CV if any([ n0<=v<=n1 for v in cell]) 
		and any([ m0<=v<=m1 for v in cell])]

# 2. Look for cells in cells12, with vertices on boundaries
cells12BV1 = [cell for cell in cells12
				if len(list(set(cell).intersection(BV1))) != 0]
cells12BV2 = [cell for cell in cells12
				if len(list(set(cell).intersection(BV2))) != 0]
pivots = sorted(AA(sorted)(cells12BV1+cells12BV2))
pivots = [cell for k,cell in enumerate(pivots[:-1]) if cell==pivots[k+1]]

# 3. Build intersection tasks
def cuttingTest(cutHyperplane,polytope):
	signs = [INNERPROD([cutHyperplane, V[v]+[1.]]) for v in polytope]
	return any([value<0 for value in signs]) and any([value>0 for value in signs])

for pivotCell in pivots:
	cutVerts = [v for v in pivotCell if v in BV]
	for v in cutVerts:
		cutFacets = VE[v]
		cells2cut = VC[v]
		print "cutFacets,cells2cut =", CART([cutFacets,cells2cut])
		for facet,cell2cut in CART([cutFacets,cells2cut]):
			polytope = CV[cell2cut]
			points = [V[w] for w in EEV[facet]]
			dim = len(points[0])
			theMat = Matrix( [(dim+1)*[1.]] + [p+[1.] for p in points] )
			cutHyperplane = [(-1)**(col)*theMat.minor(0,col).determinant() 
								for col in range(dim+1)]
			if cuttingTest(cutHyperplane,polytope):
				print "\nfacet,cell2cut =",facet,cell2cut

