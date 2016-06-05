from larlib import *
%run out_cl_alberto.lar
FV = [[v-1 for v in f] for f in FV if not all([v<=8 for v in f])]
vdict = {tuple(v):k for k,v in enumerate(V)}


def triangle2quad(vdict):
	def triangle2quad0(triangle):
		triplesOfCoords = [V[v] for v in triangle]
		(x,y,z),(X,Y,Z) = CONS([AA(min),AA(max)])(TRANS(triplesOfCoords))
		if x==X: quad = [(x,y,z),(x,Y,z),(x,Y,Z),(x,y,Z)]
		elif y==Y: quad = [(x,y,z),(X,y,z),(X,y,Z),(x,y,Z)]
		elif z==Z: quad = [(x,y,z),(X,y,z),(X,Y,z),(x,Y,z)]
		return tuple(sorted([vdict[coords] for coords in quad]))
	return triangle2quad0

def larCuboidsFacets(cells,dim=2):
   n = int(2**(dim-1))
   facets = []
   for cell in cells:
      coords = [AR([V[v],v]) for v in cell] # decorate coords with vertex index
      doubleFacets = [sorted(coords,key=(lambda x: x[k])) for k in range(dim)]
      facets += AA(AA(LAST))(CAT([[pair[:n],pair[n:]] for pair in doubleFacets]))
   return sorted(set(AA(tuple)(facets))) # remove duplicates
t1 = time.clock()
_=larCuboidsFacets(FV)
t2 = time.clock()
print "method 0 =", t2-t1

def larCuboidsFacets(FV):
	triplesOfCoords = [[V[v] for v in quad] for quad in FV]
	maxmins = AA(COMP([ CONS([ AA(min),AA(max) ]),TRANS ]))(triplesOfCoords)
	def edge2verts(args): 
		coords1,coords2 = args
		return sorted([vdict[coords1],vdict[coords2]])
	edges = set()
	for (x,y,z),(X,Y,Z) in maxmins:
		if x==X: edges = edges.union(AA(tuple)(AA(edge2verts)( 
			[((x,y,z),(x,Y,z)),((x,Y,z),(x,Y,Z)),((x,Y,Z), (x,y,Z)),((x,y,Z),(x,y,z))])))
		elif y==Y: edges = edges.union(AA(tuple)(AA(edge2verts)( 
			[((x,y,z),(X,y,z)),((X,y,z),(X,y,Z)),((X,y,Z), (x,y,Z)),((x,y,Z),(x,y,z))])))
		elif z==Z: edges = edges.union(AA(tuple)(AA(edge2verts)( 
			[((x,y,z),(X,y,z)),((X,y,z),(X,Y,z)),((X,Y,z), (x,Y,z)),((x,Y,z),(x,y,z))])))
	return sorted(edges)
t1 = time.clock()
_=larCuboidsFacets(FV)
t2 = time.clock()
print "method 1 =", t2-t1

def larCuboidsFacets(FV):
	triplesOfCoords = [[V[v] for v in quad] for quad in FV]
	maxmins = AA(COMP([ CONS([ AA(min),AA(max) ]),TRANS ]))(triplesOfCoords)
	out = set(CAT([AA(tuple)(AA(edge2verts)( 
		[((x,y,z),(x,Y,z)),((x,Y,z),(x,Y,Z)),((x,Y,Z), (x,y,Z)),((x,y,Z),(x,y,z))]))
	if x==X else AA(tuple)(AA(edge2verts)( 
		[((x,y,z),(X,y,z)),((X,y,z),(X,y,Z)),((X,y,Z), (x,y,Z)),((x,y,Z),(x,y,z))]))
	if y==Y else AA(tuple)(AA(edge2verts)( 
		[((x,y,z),(X,y,z)),((X,y,z),(X,Y,z)),((X,Y,z), (x,Y,z)),((x,Y,z),(x,y,z))]))
	for (x,y,z),(X,Y,Z) in maxmins]))
	return sorted(out)
t1 = time.clock()
_=larCuboidsFacets(FV)
t2 = time.clock()
print "method 2 =", t2-t1


FV = list(set(AA(triangle2quad(vdict))(FV)))
VIEW(STRUCT(MKPOLS((V,FV))))

EV = larCuboidsFacets(FV)
VE = invertRelation(EV)
VV = [[u for e in v for u in EV[e] if u!=k] for k,v in enumerate(VE)]


W = V
for i in range(10):
    W = [CCOMB([W[w] for w in VV[k]]) for k,v in enumerate(W)]

VIEW(STRUCT(MKPOLS((V,EV))))

csr_mat = larBoundary(FV,EV)  # 59 secs
assert len(EV)*2 == csr_mat.shape[0]*2 == csr_mat.nnz

assert AA(max)(TRANS(V)) == [261, 299, 50]

BF = set([k for k,f in enumerate(FV) if all([V[v][2]==50 for v in f])])
VIEW(SKEL_1(STRUCT(MKPOLS((V,[FV[k] for k in BF ])))))
BE = set([])

t1 = time.clock()
f_col = csr_matrix((len(FV),1),dtype='b')
for h in BF: f_col[h,0] = 1
t2 = time.clock()
print "f_col =", t2-t1

t1 = time.clock()
f_col = coo_matrix(([1]*len(BF),(list(BF),[0]*len(BF))),shape=(len(FV),1),dtype='b').tocsr()
t2 = time.clock()
print "f_col =", t2-t1
t1 = time.clock()
e_col = csr_mat * f_col
t2 = time.clock()
print "e_col =", t2-t1
edges = [e for e in e_col.nonzero()[0].tolist() if e_col[e,0]==1 ]  
delta_edges = edges
VIEW(STRUCT(MKPOLS((V,[EV[e] for e in delta_edges]))))
for i in range(50):
	# e_row = csr_matrix((1,len(EV)),dtype='b')
	# for k in delta_edges: e_row[0,k] = 1
	e_row = coo_matrix(( [1]*len(delta_edges), ([0]*len(delta_edges), delta_edges) ),
				shape=(1,len(EV)), dtype='b').tocsr()
	f_row = e_row * csr_mat 
	BF = BF.union([f for f in f_row.nonzero()[1].tolist()]) 

	# f_col = csr_matrix((len(FV),1),dtype='b')
	# for h in BF: f_col[h,0] = 1
	f_col = coo_matrix(([1]*len(BF),(list(BF),[0]*len(BF))),shape=(len(FV),1),dtype='b').tocsr()
	e_col = csr_mat * f_col
	delta_edges = [ e for e in e_col.nonzero()[0].tolist() if e_col[e,0]==1 ]
	if i%10==9: edges += delta_edges

VIEW(STRUCT(MKPOLS((W,[EV[e] for e in edges]))))
VIEW(STRUCT(MKPOLS((W,[FV[f] for f in BF]))))

TV = CAT([[(FV[f][0],FV[f][1],FV[f][2]),(FV[f][2],FV[f][1],FV[f][3])]  for f in BF])
VIEW(STRUCT(MKPOLS((V,TV))))
