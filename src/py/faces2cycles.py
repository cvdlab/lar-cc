from larlib import *
%run out_cl_alberto.lar

t1 = time.clock()
FV = [[v-1 for v in f] for f in FV]
vdict = {tuple(v):k for k,v in enumerate(V)}
t2 = time.clock()
print "\n> vertex dictionary: time =",t2-t1

#
# Transformation from triangles to (ortho-)quads
#
def triangle2quad(vdict):
	def triangle2quad0(triangle):
		triplesOfCoords = [V[v] for v in triangle]
		(x,y,z),(X,Y,Z) = CONS([AA(min),AA(max)])(TRANS(triplesOfCoords))
		if x==X: quad = [(x,y,z),(x,Y,z),(x,Y,Z),(x,y,Z)]
		elif y==Y: quad = [(x,y,z),(X,y,z),(X,y,Z),(x,y,Z)]
		elif z==Z: quad = [(x,y,z),(X,y,z),(X,Y,z),(x,Y,z)]
		return tuple(sorted([vdict[coords] for coords in quad]))
	return triangle2quad0

if __name__=="__main__":
	t1 = time.clock()
	FV = list(set(AA(triangle2quad(vdict))(FV)))
	t2 = time.clock()
	print "\n> transformation to quads: time =",t2-t1

VIEW(STRUCT(MKPOLS((V,FV))))

#
# Computation of edges of quads
#
def larQuadEdges(cells,dim=3):
   n,out = int(2**(dim-1)),[]
   for cell in cells:
      facets = []
      coords = [AR([V[v],v]) for v in cell] # decorate coords with vertex index
      doubleFacets = [sorted(coords,key=(lambda x: x[k])) for k in range(dim)]
      facets += AA(AA(LAST))(CAT([[pair[:n],pair[n:]] for pair in doubleFacets]))
      out += CAT([[tuple(face[:n/2]),tuple(face[n/2:])] for face in facets if face!=[]])
   return sorted(set(out)) # remove duplicates

if __name__=="__main__":
	t1 = time.clock()
	EV = larQuadEdges(FV)
	VE = invertRelation(EV)
	VV = [[u for e in edges for u in EV[e] if u!=k] for k,edges in enumerate(VE)]
	t2 = time.clock()
	print "\n> Computation of topological relations: time =",t2-t1

	VIEW(STRUCT(MKPOLS((V,EV))))

#
# Optimized boundary matrix computation 
#
def larBoundary(FV,EV):
	m,n = len(FV),len(EV)
	cooEF = (csrCreate(EV)*csrCreate(FV).T).tocoo()
	data = 2*[1]*n
	row = CAT([2*[k] for k in range(n)])
	col = [cooEF.col[k] for k in range(len(cooEF.data)) if cooEF.data[k]==2]
	return coo_matrix((data,(row,col)),shape=(n,m),dtype='b').tocsr()

if __name__=="__main__":
	t1 = time.clock()
	csr_mat = larBoundary(FV,EV)  
	t2 = time.clock()
	print "\n> Computation of boundary matrix: time =",t2-t1
	assert len(EV)*2 == csr_mat.shape[0]*2 == csr_mat.shape[1]*4 == csr_mat.nnz

#
# Example of brick boundary discovery 
#
max_x,max_y,max_z = AA(max)(TRANS(V))
min_x,min_y,min_z = AA(min)(TRANS(V))

#
# Initialisation of boundary cycle BF (on brick face) 
#
BF = set([k for k,f in enumerate(FV) if 
	all([V[v][0]==max_x for v in f]) or
	all([V[v][1]==max_y for v in f]) or
	all([V[v][2]==max_z for v in f]) or
	all([V[v][0]==min_x for v in f]) or
	all([V[v][1]==min_y for v in f]) or
	all([V[v][2]==min_z for v in f])
	])
VIEW(SKEL_1(STRUCT(MKPOLS((V,[FV[k] for k in BF ])))))

#
# Iterative computation of the set of closed cycles on the boundary surface 
#
def boundaryCycles(BF,FV,EV,csr_mat):
	f_col = coo_matrix(([1]*len(BF),(list(BF),[0]*len(BF))),shape=(len(FV),1),dtype='b').tocsr()
	e_col = csr_mat * f_col
	edgeStrata = [[e for e in e_col.nonzero()[0].tolist() if e_col[e,0]==1 ]]
	delta_edges = edgeStrata[0]
	for i in range(100):
		e_row = coo_matrix(( [1]*len(delta_edges), ([0]*len(delta_edges), delta_edges) ),
					shape=(1,len(EV)), dtype='b').tocsr()
		f_row = e_row * csr_mat 
		BF = BF.union([f for f in f_row.nonzero()[1].tolist()]) 

		f_col = coo_matrix(([1]*len(BF),(list(BF),[0]*len(BF))),shape=(len(FV),1),dtype='b').tocsr()
		e_col = csr_mat * f_col
		delta_edges = [ e for e in e_col.nonzero()[0].tolist() if e_col[e,0]==1 ]
		if i%10 == 9: edgeStrata += [delta_edges]
	return edgeStrata,BF

if __name__=="__main__":
	t1 = time.clock()
	edgeStrata,BF = boundaryCycles(BF,FV,EV,csr_mat)
	t2 = time.clock()
	print "\n> Computation of boundary cycles: time =",t2-t1

#
# Iterative computation of smoothed vertex positions 
#
def LaplacianSmoothing(V,VV,nsteps=10):
	W = V
	for i in range(nsteps):
		W = [CCOMB([W[w] for w in VV[k]]) for k,v in enumerate(W)]
	return W

if __name__=="__main__":
	nsteps = 10
	t1 = time.clock()
	W = LaplacianSmoothing(V,VV,nsteps)
	t2 = time.clock()
	print "\n> Laplacian smoothig: nsteps=",nsteps, ", time =",t2-t1

#
# Visualisation of set of cycles and partial surfaces 
#
if __name__=="__main__":
	VIEW(STRUCT(MKPOLS((V,[EV[e] for e in CAT(edgeStrata)]))))
	
	colors = [CYAN,MAGENTA,WHITE,RED,YELLOW,GREEN,GRAY,ORANGE,BLACK,BLUE,PURPLE,BROWN]
	VIEW(STRUCT([COLOR(colors[k%12])(STRUCT(MKPOLS((V,[EV[e] for e in curve])))) 
		for k,curve in enumerate(edgeStrata)]))
		
	layers = [COLOR(colors[k%12])(STRUCT(MKPOLS((V,[EV[e] for e in curve])))) 
		for k,curve in enumerate(edgeStrata)]
			
	for h in range(1,len(layers)): VIEW(STRUCT(layers[:h]))
		
		
	VIEW(STRUCT(MKPOLS((W,[FV[f] for f in BF]))))

#
# Visualisation of triangulated partial surfaces 
#
if __name__=="__main__":
	TV = CAT([[(FV[f][0],FV[f][1],FV[f][2]),(FV[f][2],FV[f][1],FV[f][3])]  for f in BF])
	VIEW(STRUCT(MKPOLS((W,TV))))

#
# Simplification of stratified surface
#

vcycleStrata = [makeCycles((V,[EV[e] for e in edgeStratum]))[0] for edgeStratum in edgeStrata]
simpleVertexCycles = [[[v for h,v in enumerate(vcycle) if h%10==0]+[vcycle[0]] 
	for k,vcycle in enumerate(vcycleStratum) if k%2==0] 
		for vcycleStratum in vcycleStrata]

simpleVertexCycles = [[cycle for cycle in level if cycle!=[1,1]] for level  in simpleVertexCycles]

VIEW(STRUCT([ COLOR(colors[k%12])( STRUCT([ POLYLINE([V[v] for v in polyline ]) 
	for polyline in level ])) for k,level in enumerate(simpleVertexCycles)]))

oldEdges = [CAT([[ (u,v) for u,v in  zip(polyline[:-1],polyline[1:])+[(polyline[-1],polyline[0])]] 
	for polyline in level if len(polyline)>2]) for level in simpleVertexCycles]

vertsPerLevel = AA(CAT)(simpleVertexCycles)
vertsPerLevel = AA(sorted)((AA(set)(vertsPerLevel)))
trees = [spatial.KDTree(AA(np.array)([V[v] for v in vlevel]))
	for vlevel in vertsPerLevel]

stripes = zip(trees[:-1],trees[1:])
newEdges = []
for k,(tree1,tree2) in enumerate(stripes):

	nearest = tree1.query(tree2.data)[1]
	verts1 = vertsPerLevel[k+1]
	verts2 = [vertsPerLevel[k][h] for h in nearest]
	newEdges += [AA(sorted)(zip(verts1,verts2))]

	nearest = tree2.query(tree1.data)[1]
	verts1 = vertsPerLevel[k]
	verts2 = [vertsPerLevel[k+1][h] for h in nearest]
	newEdges += [AA(sorted)(zip(verts1,verts2))]
	
newEdges = sorted(set(CAT(newEdges)))

VIEW(STRUCT(AA(POLYLINE)([[V[u],V[v]] for u,v in newEdges+CAT(oldEdges)])))



