from larlib import *
%run out_cl_alberto.lar
t1 = time.clock()
FV = [[v-1 for v in f] for f in FV if not all([v<=8 for v in f])]
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
assert AA(max)(TRANS(V)) == [261, 299, 50]

#
# Initialisation of boundary cycle BF (on brick face) 
#
BF = set([k for k,f in enumerate(FV) if all([V[v][2]==50 for v in f])])
VIEW(SKEL_1(STRUCT(MKPOLS((V,[FV[k] for k in BF ])))))

#
# Iterative computation of the set of closed cycles on the boundary surface 
#
def boundaryCycles(BF,FV,EV,csr_mat):
	f_col = coo_matrix(([1]*len(BF),(list(BF),[0]*len(BF))),shape=(len(FV),1),dtype='b').tocsr()
	e_col = csr_mat * f_col
	edges = [e for e in e_col.nonzero()[0].tolist() if e_col[e,0]==1 ]  
	delta_edges = edges
	for i in range(50):
		e_row = coo_matrix(( [1]*len(delta_edges), ([0]*len(delta_edges), delta_edges) ),
					shape=(1,len(EV)), dtype='b').tocsr()
		f_row = e_row * csr_mat 
		BF = BF.union([f for f in f_row.nonzero()[1].tolist()]) 

		f_col = coo_matrix(([1]*len(BF),(list(BF),[0]*len(BF))),shape=(len(FV),1),dtype='b').tocsr()
		e_col = csr_mat * f_col
		delta_edges = [ e for e in e_col.nonzero()[0].tolist() if e_col[e,0]==1 ]
		if i%10 == 9: edges += delta_edges
	return edges,BF

if __name__=="__main__":
	t1 = time.clock()
	edges,BF = boundaryCycles(BF,FV,EV,csr_mat)
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
	VIEW(STRUCT(MKPOLS((W,[EV[e] for e in edges]))))
	VIEW(STRUCT(MKPOLS((W,[FV[f] for f in BF]))))

#
# Visualisation of triangulated partial surfaces 
#
if __name__=="__main__":
	TV = CAT([[(FV[f][0],FV[f][1],FV[f][2]),(FV[f][2],FV[f][1],FV[f][3])]  for f in BF])
	VIEW(STRUCT(MKPOLS((V,TV))))

