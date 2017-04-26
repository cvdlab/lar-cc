from larlib import *

#
# personalised camera
#

def MYVIEW(obj):
  glcanvas=GLCanvas()
  glcanvas.setOctree(Octree(Plasm.getBatches(obj)))
  glcanvas.frustum.pos=Vec3f(3.482283e+02,9.733051e+01,1.291512e+02)
  glcanvas.frustum.dir=Vec3f(-7.539707e-01,2.843412e-01,-5.921810e-01)
  glcanvas.frustum.vup=Vec3f(-5.159718e-01,3.015976e-01,8.017555e-01)
  glcanvas.redisplay()
  glcanvas.runLoop()

def MYVIEW(obj):
  glcanvas=GLCanvas()
  glcanvas.setOctree(Octree(Plasm.getBatches(obj)))
  glcanvas.frustum.pos=Vec3f(2.406189e+02,8.040689e+01,1.165992e+02)
  glcanvas.frustum.dir=Vec3f(-5.930701e-02,4.281818e-01,-9.017444e-01)
  glcanvas.frustum.vup=Vec3f(-2.249289e-01,8.743718e-01,4.299777e-01)
  glcanvas.redisplay()
  glcanvas.runLoop()

def MYVIEW(obj):
  glcanvas=GLCanvas()
  glcanvas.setOctree(Octree(Plasm.getBatches(obj)))
  glcanvas.frustum.pos=Vec3f(1.432074e+02,-3.626540e+01,1.961070e+02)
  glcanvas.frustum.dir=Vec3f(-4.973759e-02,7.350764e-01,-6.761576e-01)
  glcanvas.frustum.vup=Vec3f(5.976993e-02,6.779727e-01,7.326531e-01)
  glcanvas.redisplay()
  glcanvas.runLoop()

def MYVIEW(obj):
  glcanvas=GLCanvas()
  glcanvas.setOctree(Octree(Plasm.getBatches(obj)))
  glcanvas.frustum.pos=Vec3f(1.743804e+02,2.692199e+01,1.217709e+02)
  glcanvas.frustum.dir=Vec3f(-4.973787e-02,7.350735e-01,-6.761605e-01)
  glcanvas.frustum.vup=Vec3f(5.976970e-02,6.779758e-01,7.326503e-01)
  glcanvas.redisplay()
  glcanvas.runLoop()

#
# data input
#
#%run 16-06-13/out_cl.lar
%run out_cl_alberto.lar

t1 = time.clock()
FV = [sorted([v-1 for v in f]) for f in FV if not any([v<8 for v in f])]
V = [v for k,v in enumerate(V) if k>7]
FV = [[v-8 for v in f] for f in FV]

vdict = {tuple(v):k for k,v in enumerate(V)}
t2 = time.clock()
print "\n> vertex dictionary: time =",t2-t1

#
# Transformation from triangles to (ortho-)quads 
#
if len(FV[0]) == 3:
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
else: pass

MYVIEW(STRUCT(MKPOLS((V,FV))))

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

	#MYVIEW(STRUCT(MKPOLS((V,EV))))

#
# Optimized boundary matrix computation 
#
def larBoundary(FV,EV):
	m,n = len(FV),len(EV)
	edict = {tuple(e):k for k,e in enumerate(EV)}
	data = 2*[1]*n
	row = []; [row.extend([ edict[(f[i],f[j])] for i in range(4) for j in range(i+1,4) 
			if (f[i],f[j]) in edict ]) for f in FV]
	col = []; [col.extend([k,k,k,k]) for k in range(m)]
	return coo_matrix((data,(row,col)),shape=(n,m),dtype='b').tocsc()

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
"""
BF = sorted([k for k,f in enumerate(FV) if 
	all([V[v][2]==max_z for v in f])
	])
BF = sorted([k for k,f in enumerate(FV) if 
	all([V[v][0]==max_x for v in f]) or
	all([V[v][1]==max_y for v in f]) or
	all([V[v][2]==max_z for v in f]) or
	all([V[v][0]==min_x for v in f]) or
	all([V[v][1]==min_y for v in f]) or
	all([V[v][2]==min_z for v in f])
	])
"""
BF = sorted([k for k,f in enumerate(FV) if 
	all([V[v][0]==max_x for v in f]) or
	all([V[v][1]==max_y for v in f]) or
	all([V[v][2]==max_z for v in f]) 
	])
	
MYVIEW(SKEL_1(STRUCT(MKPOLS((V,[FV[k] for k in BF ])))))

#
# Iterative computation of the set of closed cycles on the boundary surface 
#
def boundaryCycles(V,BF,FV,EV,csr_mat):
	V = AA(tuple)(V)
	lenEV,lenFV = len(EV),len(FV)
	f_col = coo_matrix(([1]*len(BF), (BF, [0]*len(BF))),
							shape=(len(FV),1),dtype='b').tocsc()
	e_col = csr_mat * f_col
	edgeStrata = [[e for e in e_col.nonzero()[0].tolist() if e_col[e,0]==1 ]]
	delta_edges = edgeStrata[0]
	delta_faces = BF
	previous_ncycles = 0
	
	MYVIEW(STRUCT(MKPOLS((V,[EV[e] for e in delta_edges]))))
	
	for i in range(110):
		if delta_edges != []:
			t1 = time.clock()
			e_len = len(delta_edges)
			e_row = coo_matrix(( [1]*e_len, ([0]*e_len, delta_edges) ),
						shape=(1,lenEV), dtype='b').tocsr()
			f_row = e_row * csr_mat 
			delta_faces = list(set(f_row.nonzero()[1]).difference(delta_faces))
			#MYVIEW(STRUCT(MKPOLS((V,[FV[f] for f in delta_faces]))))
		
			f_col = coo_matrix(([1]*len(delta_faces), (delta_faces, 
					[0]*len(delta_faces))),	shape=(len(FV),1),dtype='b').tocsc()
			e_col = csr_mat * f_col
			newEdges = [e for e in e_col.nonzero()[0] if e_col[e,0]==1 ]
			#MYVIEW(STRUCT(MKPOLS((V,[EV[e] for e in newEdges]))))
			delta_edges = list(set(newEdges).difference(delta_edges))
			"""
			G=nx.Graph()
			G.add_nodes_from(V)
			G.add_edges_from([EV[e] for e in delta_edges])
			"""
			#vcycles = makeCycles((V,[EV[e] for e in delta_edges]))[0] 
			#ncycles = nx.number_connected_components(G)
			#print "> ncycles =",ncycles
		
			t2 = time.clock()
			print "> i =",i, ", time =",t2-t1		
		
			if i%11 == 10: #or ncycles != previous_ncycles: 
				if delta_edges!=[]: edgeStrata += [delta_edges]
				#previous_ncycles = ncycles

	MYVIEW(STRUCT(MKPOLS((V,[EV[e] for e in CAT(edgeStrata)]))))
	return edgeStrata,BF

if __name__=="__main__":
	t1 = time.clock()
	edgeStrata,BF = boundaryCycles(V,BF,FV,EV,csr_mat)
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
	MYVIEW(STRUCT(MKPOLS((V,[EV[e] for e in CAT(edgeStrata)]))))
	
	colors = [CYAN,MAGENTA,WHITE,RED,YELLOW,GREEN,GRAY,ORANGE,BLACK,BLUE,PURPLE,BROWN]
	MYVIEW(STRUCT([COLOR(colors[k%12])(STRUCT(MKPOLS((V,[EV[e] for e in curve])))) 
		for k,curve in enumerate(edgeStrata)]))
		
	layers = [COLOR(colors[k%12])(STRUCT(MKPOLS((V,[EV[e] for e in curve])))) 
		for k,curve in enumerate(edgeStrata)]
			
	for h in range(1,len(layers)): MYVIEW(STRUCT(layers[:h]))
		
		
	MYVIEW(STRUCT(MKPOLS((W,[FV[f] for f in BF]))))

#
# Visualisation of triangulated partial surfaces 
#
if __name__=="__main__":
	TV = CAT([[(FV[f][0],FV[f][1],FV[f][2]),(FV[f][2],FV[f][1],FV[f][3])]  for f in BF])
	MYVIEW(STRUCT(MKPOLS((W,TV))))

#
# Simplification of stratified surface
#

vcycleStrata = [makeCycles((V,[EV[e] for e in edgeStratum]))[0] for edgeStratum in edgeStrata]
simpleVertexCycles = [[[v for h,v in enumerate(vcycle) if h%10==0]+[vcycle[0]] 
	for k,vcycle in enumerate(vcycleStratum) if k%2==0] 
		for vcycleStratum in vcycleStrata]

simpleVertexCycles = [[cycle for cycle in level if cycle!=[1,1]] for level  in simpleVertexCycles]

MYVIEW(STRUCT([ COLOR(colors[k%12])( STRUCT([ POLYLINE([V[v] for v in polyline ]) 
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
	""" 
	nearest = tree2.query(tree1.data)[1]
	verts1 = vertsPerLevel[k]
	verts2 = [vertsPerLevel[k+1][h] for h in nearest]
	newEdges += [AA(sorted)(zip(verts1,verts2))]
	""" 
newEdges = sorted(set(AA(tuple)(CAT(newEdges))))

MYVIEW(STRUCT(AA(POLYLINE)([[V[u],V[v]] for u,v in newEdges])))
MYVIEW(STRUCT(AA(POLYLINE)([[V[u],V[v]] for u,v in newEdges+CAT(oldEdges)])))

# with Laplacian smoothing
MYVIEW(STRUCT(AA(POLYLINE)([[W[u],W[v]] for u,v in newEdges+CAT(oldEdges)])))



