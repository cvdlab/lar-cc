from larlib import *
from meshpy.tet import MeshInfo, build, Options

# LAR model with non-contractible faces
# ------------------------------------------------------------------------------

V = [[0.25, 0.25, 0.0], [0.25, 0.75, 0.0], [0.75, 0.75, 0.0], [0.75, 0.25, 0.0], [1.0, 
0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.25, 0.25, 1.0], [0.25, 
0.25, 2.0], [0.25, 0.75, 2.0], [0.25, 0.75, 1.0], [0.25, 0.75, -1.0], [0.25, 0.25, 
-1.0], [0.75, 0.75, -1.0], [0.75, 0.25, -1.0], [0.75, 0.25, 1.0], [0.75, 0.75, 1.0], 
[1.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0], [0.75, 0.75, 2.0], 
[0.75, 0.25, 2.0]]

CV = [(0,1,2,3,4,5,6,7,8,11,16,17,18,19,20,21), 
(0,1,2,3,8,11,16,17),
(0,1,2,3,12,13,14,15), 
(8,9,10,11,16,17,22,23)]

FV = [(2,3,16,17),(6,7,20,21),(12,13,14,15),(0,1,8,11),(1,2,11,17),(0,1,12,13),
(4,6,18,20),(5,7,19,21),(0,3,13,15),(0,3,8,16),(0,1,2,3),
(10,11,17,22),(2,3,14,15),(8,9,16,23),(8,11,16,17),
(1,2,12,14),(16,17,22,23),(4,5,18,19),(8,9,10,11),(
9,10,22,23),(0,1,2,3,4,5,6,7),(8, 11,16,17,18,19,20,21)]

EV =[(3,15),(7,21),(10,11),(4,18),(12,13),(5,19),(8,9),(18,19),(22,23),(0,3),(1,11),
(16,17),(0,8),(6,7),(20,21),(3,16),(10,22),(18,20),(19,21),(1,2),(12,14),(4,5),(
8,11),(13,15),(16,23),(14,15),(11,17),(17,22),(2,14),(2,17),(0,1),(9,10),(8,16),
(4,6),(1,12),(5,7),(0,13),( 9,23),(6,20),(2,3)]

VV = AA(LIST)(range(len(V)))
hpc = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV,CV],hpc,0.6))

BF = boundaryCells(CV,FV)
VIEW(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,[FV[f] for f in BF],EV))))

VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV+EV+VV))))
VIEW(EXPLODE(1.2,1.2,1.2)(AA(SKEL_1)(MKPOLS((V,CV)))))


# Correction of non-signed boundary op for general (non contractible) LAR cells 
# ------------------------------------------------------------------------------

def boundary3(CV,FV):
    lenV = max(CAT(CV))+1
    FC = invertRelation(crossRelation(lenV,CV,FV))
    oddAdjacencies = [f for f,cells in enumerate(FC) if len(cells)>2]
    if oddAdjacencies == []: 
        return boundary(CV,FV)
    else:
        FF = crossRelation(lenV,FV,FV)
        for f in oddAdjacencies:
            anomalies = list(set(FF[f]).difference([f]))
            for g in anomalies:
                if {f}.issubset(FF[g]): 
                    FC[g]= list(set(FC[f]).difference(FC[g]))
                    FC[f]= list(set(FC[f]).difference(FC[g]))
        out = csr_matrix((len(FV),len(CV)),dtype='b')
        for h in range(len(FV)):
            for k in FC[h]:
                out[h,k] = 1
        return out

V,[VV,EV,FV,CV] = larCuboids([2,1,1],True)
mod1 = Struct([(V,FV,EV),t(.25,.25,0),s(.25,.5,2),(V,FV,EV)])
V,FV,EV = struct2lar(mod1)

W,FW,EW = V,FV,EV
quadArray = [[W[v] for v in face] for face in FW]
parts = boxBuckets3d(containmentBoxes(quadArray))
Z,FZ,EZ = spacePartition(W,FW,EW, parts)
Z,FZ,EZ = larSimplify((Z,FZ,EZ),radius=0.0001)
V,FV,EV = Z,FZ,EZ

CF = AA(sorted)([[20,12,21,5,19,6],
[27,1,5,28,13,23],
[12,14,25,17,10,4],
[1,7,17,24,11,18],
[30,29,26,16,8,22,10,11,4,18,24,25],
[2,3,8,9,0,15]])

CV = [list(set(CAT([FV[f]  for f in faces]))) for faces in CF]

CV = [[10, 11, 12, 13, 18, 19, 20, 21],
 [18, 19, 20, 21, 22, 23, 25, 26],
 [0, 1, 4, 5, 10, 13, 18, 21],
 [2, 3, 4, 5, 18, 21, 25, 26],
 [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 16, 17, 18, 21, 24, 25, 26, 27],
 [6, 8, 14, 15, 16, 24, 28, 29]]

VV = AA(LIST)(range(len(V)))
hpc = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV,CV],hpc,0.6))

BF = boundaryCells(CV,FV)
VIEW(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,[FV[f] for f in BF],EV))))


# Tetrahedralization 
# ------------------------------------------------------------------------------

CF = crossRelation(len(V),CV,FV)

FE = crossRelation(len(V),FV,EV)  # correct (for general LAR cells)
cycles = []
for faceEdges in FE:
    vcycles,_ = makeCycles((V,[EV[e] for e in faceEdges]))
    cycles += [[vcycle for k,vcycle in enumerate(vcycles) if k%2==0]]
    
cycles = [[[16, 17, 2, 3]],   # removed dups ...
 [[20, 21, 7, 6]],
 [[14, 15, 13, 12]],
 [[8, 11, 1, 0]],
 [[11, 17, 2, 1]],
 [[12, 13, 0, 1]],
 [[18, 20, 6, 4]],
 [[19, 21, 7, 5]],
 [[13, 15, 3, 0]],
 [[8, 16, 3, 0]],
 [[17, 22, 10, 11]],
 [[14, 15, 3, 2]],
 [[16, 23, 9, 8]],
 [[12, 14, 2, 1]],
 [[22, 23, 16, 17]],
 [[18, 19, 5, 4]],
 [[10, 11, 8, 9]],
 [[22, 23, 9, 10]],
 [[6, 7, 5, 4], [2, 3, 0, 1]],
 [[20, 21, 19, 18], [16, 17, 11, 8]]]



CF = [[1, 7, 17, 6, 0, 4, 9, 3, 20,21], # plus 20,21
 [14, 0, 4, 10, 9, 3],  # no 20,21
 [2, 12, 15, 10, 8, 5],        # no 20
 [16, 11, 19, 18, 14, 13]]   # no 21
 
 
    
VIEW(EXPLODE(1.2,1.2,1.2)(AA(POLYLINE)([[V[v] for v in cycle]+[V[cycle[0]]] 
     for cycle in CAT(cycles)])))

def faces(tet):
	v1,v2,v3,v4 = tet
	return [(v2,v3,v4),(v3,v1,v4),(v1,v2,v4),(v2,v1,v3)]

def edges(tria):
	v1,v2,v3 = tria
	return AA(sorted)([(v1,v2),(v1,v3),(v2,v3)])


def brep2lar(model,cycles,holes):
    V,FV,EV = model    
    FE = crossRelation(V,FV,EV)
    mesh_info = MeshInfo()
    mesh_info.set_points(V)
    mesh_info.set_facets_ex(cycles)
    #mesh_info.set_holes(holes)
    mesh = build(mesh_info,options=Options("pqYY"))
    W = [v for h,v in enumerate(mesh.points)]
    CW = [tet for k,tet in enumerate(mesh.elements)]
    
    def simplify(fun):
        def simplify0(simplices):
	    cellDict = defaultdict(tuple)
	    for cell in CAT(AA(fun)(simplices)):
	        cellDict[tuple(sorted(cell))] = tuple(cell)
	    return cellDict.values()
	return simplify0
    
    FW = sorted(simplify(faces)(CW))
    EW = sorted(simplify(edges)(FW))
    return W,CW,FW,EW

holes = [[0.5,0.5,0.5]]
W,CW,FW,EW = brep2lar((V,FV,EV),cycles,holes)  
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,FW))))

BF = signedSimplicialBoundary(CW,FW)
bf = (BF * mat(len(CW)*[1]).T).tolist()
bfs = [k for k,face in enumerate(CAT(bf)) if ABS(face)==1]
VIEW(EXPLODE(2,2,2)(MKPOLS((W,[FW[f] for f in bfs]))))



# Tetrahedralization 
# ------------------------------------------------------------------------------


frame = boundary(CV,FV)
CVW = crossRelation(len(V),CV,CW)  # tetrahedra by LAR cell  # maps cells to tetrahedra
CVW[0] = list(set(CVW[0]).difference(CVW[1]))  # removed double tetrahedra

tetraBreps = [[faces(CW[t]) for t in c] for c in CVW]
[Volume([W,tetra]) for tetra in CAT(tetraBreps)]  #  test of tetrahedra coherent orientation 
for c in range(len(CV)):
	VIEW(STRUCT(MKPOLS((V,CAT(tetraBreps[c])))))   # 2-scheletro di ogni cella
    
BCt = []
for c in range(len(CV)):
	A = [CW[t] for t in CVW[c]]
	B = CAT(tetraBreps[c])
	cellBoundaryChain = boundary(A,B)*mat(len(A)*[[1]])
	bc = [k for k,val in enumerate (cellBoundaryChain) if val==1]
	print "c,bc =",c,bc
	BCt += [bc]
	VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,[B[t] for t in bc]))))  # 2-bordo di ogni cella


FVt = [[t for t in FW if set(t).issubset(f)] for f in FV]
FVt
[[(16, 2, 3), (16, 17, 2)],
 [(6, 7, 21), (21, 20, 6)],
 [(12, 15, 13), (15, 12, 14)],
 [(0, 1, 11), (11, 8, 0)],
 [(1, 17, 2), (1, 17, 11)],
 [(0, 1, 12), (0, 12, 13)],
 [(4, 20, 18), (20, 4, 6)],
 [(5, 21, 7), (21, 5, 19)],
 [(15, 0, 13), (15, 3, 0)],
 [(0, 16, 3), (16, 0, 8)],
 [(0, 1, 2), (3, 2, 0)],
 [(10, 17, 11), (17, 10, 22)],
 [(2, 15, 14), (15, 2, 3)],
 [(8, 23, 9), (23, 8, 16)],
 [(17, 8, 11), (17, 8, 16)],
 [(1, 14, 12), (14, 1, 2)],
 [(17, 22, 16), (23, 16, 22)],
 [(4, 19, 5), (19, 4, 18)],
 [(8, 10, 11), (9, 10, 8)],
 [(10, 9, 22), (23, 22, 9)],
[(0, 1, 2), (0, 5, 7), (0, 7, 1), (1, 6, 2), (3, 2, 0), (3, 2, 4), (4, 0, 3), (4, 2, 6), (4, 5, 0), (6, 1, 7)],
[(8, 21, 19), (11, 20, 21), (16, 18, 20), (17, 8, 11), (17, 8, 16), (17, 16, 20), (17, 20, 11), (18, 8, 19), (18, 16, 8), (21, 8, 11)]]

    

ET = crossRelation(len(W),EV,FW)
EVW = extendEV(EV,ET,FW)
triaModel, larModel = (FW,EVW), (FV,EV)
FT = crossRelation(len(V),FV,FW)
boundary2op = larSignedBoundary(larModel,triaModel,FT)
boundary2op.todense()

m,n,p = AA(len)([CV,FV,EV])
absFE = []
for f,face in enumerate(FV):
    faceVect = zeros((n,1))
    faceVect[f] = [1]
    edgeVect = CAT((boundary2op * faceVect).tolist())
    absFE += [[int(value)*e for e,value in enumerate(edgeVect) if value!=0]]

absFE
[[11, 15, -29, 39],
 [1, 13, -14, -38],
 [-4, 20, -23, 25],
 [10, -12, -22, 30],
 [-10, -19, -26, -29],
 [4, 30, 34, -36],
 [-3, -17, 33, 38],
 [-1, 5, 18, -35],
 [0, -9, 23, 36],
 [-9, 12, -15, 32],
 [9, 19, 30, -39],
 [-2, 16, -26, -27],
 [0, -25, -28, 39],
 [-6, 24, 32, -37],
 [11, 22, 26, 32],
 [19, -20, 28, -34],
 [8, 11, -24, 27],
 [3, -5, 7, -21],
 [2, 6, -22, 31],
 [-8, -16, -31, 37],
 [9, -13, 19, 21, 30, -33, 35, -39],
 [-7, 11, 14, 17, -18, 22, 26, 32]]


absCt = []
for c,cell in enumerate(CVW):
    cellVect = zeros((len(CW),1))
    for T in CVW[c]: cellVect[T] = [1]
    trianglesVect = CAT((BF * cellVect).tolist())
    absCt += [[int(value)*t for t,value in enumerate(trianglesVect) if abs(value)==1]]

for k in range(len(absCt)):
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,[FW[abs(t)] for t in absCt[k]]))))


CF = crossRelation(len(V), CV,FV)

[k for k,val in enumerate(CAT((boundary(CV,FV)*(mat([[0,1,0,0]])).T).tolist())) if val!=0]


for f in range(frame.shape[0]):
    row = frame[f]
    for c in range(frame.shape[1]):
        if frame[f,c] != 0:
            print f,c, frame[f,c]
            
