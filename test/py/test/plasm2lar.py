from larlib import *


def curation(strange,Z,CV,FV,EV):
	v_,ee_ = TRANS(strange)
	e_ = set(CAT(ee_))
	vv = set(CAT([EV[e] for e in e_]))
	VIEW(STRUCT(AA(POLYLINE)([[EV[e1],EV[e2]] for e1,e2 in ee_])))
	import networkx as nx
	G=nx.Graph()
	G.add_nodes_from(e_)
	G.add_edges_from(ee_)
	nx.connected_components(G)
	edges2remove = []; newEdges = []
	for edges in nx.connected_components(G):
		print edges
		edict = defaultdict(list)
		verts = CAT([EV[e] for e in edges])
		for v in verts: 
			edict[v] += [1]
		newEdge = [v for v in edict.keys() if len(edict[v])==1]
		edges2remove += list(edges)
		newEdges += [newEdge]
	v = 0; 
	W = [None for k in range(len(Z)-len(v_))]
	newindex = [None for k in range(len(Z))]
	for k,vert in enumerate(Z):
		if k not in v_:
			W[v] =  vert
			newindex[k]=v
			v += 1
	CW = [list(set([newindex[v] for v in cv if newindex[v]!=None])) for cv in CV]
	FW = [[newindex[v] for v in fv if newindex[v]!=None] for fv in FV if len(set(fv))>1]
	print "FV = ",FV
	EW = [[newindex[v] for v in ev if newindex[v]!=None] for e,ev in enumerate(EV) if e not in edges2remove]
	EW = [[w1,w2] for w1,w2 in EW if w1!=w2]  + [[newindex[v1],newindex[v2]] for [v1,v2] in newEdges]
	return W,CW,FW,EW

def hpc2lar(hpc):
	V, cells, _ = UKPOL(hpc)
	W, faces, _ = UKPOL(SKEL_2(hpc))
	U, edges, _ = UKPOL(SKEL_1(hpc))
	vertdict = defaultdict(list)
	for k,v in enumerate(V+W+U):
		vertdict[vcode(4)(v)] = k+1
	vert = dict(zip(vertdict.keys(),range(len(vertdict))))
	Z = [eval(p) for p,k in sorted(vert.items(),key=S2)]
	CV = list(set([tuple(sorted([vert[vcode(4)(V[v-1])] for v in cell])) for cell in cells]))
	FV = list(set([tuple(sorted([vert[vcode(4)(W[v-1])] for v in face])) for face in faces]))
	EV = list(set([tuple(sorted([vert[vcode(4)(U[v-1])] for v in edge])) for edge in edges]))

	VE = invertRelation(EV)
	strange = [(v,ve) for v,ve in enumerate(VE) if len(ve)==2]
	# TODO: generalize for multiple separated sequences in "strange"
	if strange != []:
		W,CW,FW,EW = curation(strange,Z,CV,FV,EV)
		return W,CW,FW,EW
	else: return Z,CV,FV,EV	

def computeBoundary(W,CW,FW):
	csrBoundaryMat = larBoundary(CW,FW)
	coo = coo_matrix(csrBoundaryMat)
	data,row,col = coo.data,coo.row,coo.col
	# boundary faces and incident cells
	fc = [(r,c) for k,(r,c) in enumerate(zip(row,col)) if csrBoundaryMat[r].nnz==1]
	boundaryFaces = []
	for face,cell in fc:
		boundaryFace = []
		h = COVECTOR([W[w] for w in FW[face]])
		for v in CW[cell]:
			if abs(INNERPROD([h,[1]+W[v]])) < 10**-7:
				boundaryFace += [v]
		boundaryFaces += [boundaryFace]
	return boundaryFaces
	
def internalFaces(W,CW):
	def brc2Csr(a): return coo2Csr(brc2Coo(a))
	CVC = brc2Csr(CW) * brc2Csr(CW).T
	# indices of 3-cells adjacent to 3-cells
	FV = []
	n = CVC.shape[1]
	for i in range(n-1):
		FVi = []
		for j in range(i+1,n):
			if CVC[i,j]>=3:
				FVi += [j]
		FV += [FVi]
	FZ = []
	for h in range(len(FV)):
		Ch = set(CW[h])
		for k in FV[h]:
			f_candidate = list(Ch.intersection(CW[k]))
			print "f_candidate",f_candidate
			verts = array([W[v] for v in f_candidate])
			centroid = CCOMB(verts)
			aligned = numpy.linalg.norm(cross(verts[1]-verts[0],-(verts[2]-centroid))) < 10**-5
			if not aligned: FZ += [f_candidate]
	return FZ

def permutationOrbits(d):
    out = []
    while d:
        x = list(d)[0]
        orbit = []
        while x in d:
            orbit += [x],
            x = d.pop(x)
        out += [CAT(orbit)+orbit[0]]
    return out

def SBoundary3(W,EW,FW):
   SB_2 = SBoundary2(EW,FW)
   D,I,J = [],[],[]
   m,n = SB_2.shape
   marks = [0 for k in range(n)]
   store = [0 for k in range(n)]
   # permutation subgroups of edges
   FE = crossRelation0(len(W),FW,EW)
   #FE = [list(SB_2[:,f].tocoo().row) for f in range(SB_2.shape[1])]
   EF_angle, ET,TW,FT = faceSlopeOrdering((W,FW,EW),FE)
   cellCount = -1

   while sum(marks) < 2*n:
      # choose f
      f = Choose(marks)
      print "f =",f
      
      # start 2-chain extraction from f seed
      if marks[f] == 0: c_2 = [(f,1)] 
      elif marks[f] == 1: c_2 = [(f,-store[f])]    
      Stripe_2 = coo_matrix(([],([],[])),(n,1))
   
      # computation of c_2 boundary
      C_2 = chain2coords(c_2,n).tocsc() + Stripe_2
      C_1 = SB_2 * C_2

      while C_1.nnz != 0:  
         stripe = dict()
         
         # computation of coboundaries of c_2 boundary
         dict_C_1 = coords2chainDict(C_1)
         for cell,code in dict_C_1.items():
            C1 = chain2coords([(cell,code)],m)
            C2 = C1.T * SB_2
            subgroup = list(C2.tocoo().col)
            pivot = (set(subgroup).intersection(C_2.tocoo().row)).pop()
            if code == 1: 
               adj = next(EF_angle[cell])(pivot)
            elif code == -1:
               adj = prev(EF_angle[cell])(pivot)
            if SB_2[cell,adj] == SB_2[cell,pivot] :
               stripe[adj] = -1 * C_2[pivot,0]
            else:
               stripe[adj] = 1 * C_2[pivot,0]
         Stripe_2 = chain2coords(stripe.items(),n).tocsc()
         C_2 += Stripe_2
         C_1 = SB_2 * C_2
      
      cellCount += 1
      facets = list(C_2.tocoo().row)
      coeffs = list(C_2.tocoo().data)
      cells = [cellCount] * len(coeffs)
      for k,facet in enumerate(facets): 
         marks[facet] += 1
         store[facet] += coeffs[k]
      D += coeffs
      I += facets
      J += cells
   return coo_matrix((D,(I,J)),dtype=int).tocsc()


def Boundary3(W,EW,FW,boundaryFaces):
	SB_3 = SBoundary3(W,EW,FW)
	vertDict = dict([(tuple(v),k) for k,v in enumerate(W)])
	vdict = sorted(vertDict)
	first, last = vertDict[vdict[0]], vertDict[vdict[-1]]
	WF = invertRelation(boundaryFaces)
	fmins, fmaxs = set(WF[first]), set(WF[last])
	CF = [list(SB_3[:,c].tocoo().row) for c in range(SB_3.shape[1])]
	exterior = [k for k,cell in enumerate(CF) if (fmins.intersection(cell) == fmins) 
		and (fmaxs.intersection(cell) == fmaxs)][0]
	# TODO: generalize the test for the other coordinates, to treat the (very) unlikely cases that this test doesn't work
	m,n = SB_3.shape
	coo = coo_matrix(SB_3)
	triples = []
	for (d,i,j) in zip(coo.data,coo.row,coo.col):
		if j < exterior:
			triples += [(d,i,j)]
		if j > exterior:
			triples += [(d,i,j-1)]
	data,row,col = TRANS(triples)
	return csc_matrix((data,(row,col)),dtype=int)

if __name__=="__main__":

	## Input definition
	
	hpc = UNION([CUBE(1),T(1)(1)(CUBE(1))])
	hpc = UNION([CUBE(1),T(1)(1.5)(CUBE(2))]) # KO !!
	hpc = MKPOL([
		[[0.5, 1.4901161193847656e-08,-0.2], [0.10000000894069672, 1.4901161193847656e-08,-0.2],
		[0.5, 1.4901161193847656e-08,0.0], [0.10000000894069672, 1.4901161193847656e-08,0.0],
		[0.5, 1.4901161193847656e-08,0.2], [0.10000000894069672, 1.4901161193847656e-08,0.2],
		[0.3535533547401428, 0.3535533845424652,0.0],[0.0707106813788414, 0.070710688829422,0.0],
		[0.3535533547401428, 0.3535533845424652,0.2],[0.0707106813788414, 0.070710688829422,0.2],
		[0.3535533547401428, 0.3535533845424652,0.4],[0.0707106813788414, 0.070710688829422,0.4],
		[0.0, 0.5,0.2],[0.0, 0.10000000894069672,0.2],
		[0.0, 0.5, 0.4], [0.0, 0.10000000894069672, 0.4]],
		[[1,2,3,4,7,8],[3,4,5,6,7,8,9,10],[7,8,9,10,13,14],[9,10,11,12,13,14,15,16]],
		None
    ])
	hpc = DIFFERENCE([CUBE(2),CUBE(1)])
	hpc = MKPOL([
		[[0,0,0],[1,0,0],[1,1,0],[0,1,0],
		[0,0,0.5], [1,0,0.5], [1,1,0.5], [0,1,0.5],
		[0,0,1],[1,0,1],[1,1,1],[0,1,1]],
		[[1,2,3,4,5,6,7,8],[5,6,7,8,9,10,11,12]],
		None
	])
	hpc = MKPOL([
		[[0,0,0],[1,0,0],[0,0,0.3],[1,0,0.3],[0,0,0.6],[1,0,0.6],
		[0,1,0.3],[1,1,0.3],[0,1,0.6],[1,1,0.6]],
		[[1,2,3,4,7,8],[3,4,5,6,7,8,9,10]],
		None
	])
	hpc = UNION([CUBE(1),T(1)(1)(CUBE(2))])
	hpc = UNION([CUBE(1),T(1)(0.5)(CUBE(2))]) # KO !!

	VIEW(hpc)
	VIEW(SKEL_1(hpc))
	
	## HPC -> LAR transformation
	
	W,CW,FW,EW = hpc2lar(hpc)
	VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,FW+EW))))

	## Boundary triangulation
	
	BF = [FW[f] for f in lar2boundaryFaces(CW,FW)]
	
	boundaryFaces = [face for face in computeBoundary(W,CW,FW) if face!=[]]
	VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,boundaryFaces))))

	FW = boundaryFaces + internalFaces(W,CW)
	FE = crossRelation0(len(W),FW,EW)
	triangleSet = boundaryTriangulation(W,boundaryFaces,EW,FE)
	TW,FT = triangleIndices(triangleSet,W)
	VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,TW))))

	## Cell numbering

	submodel = STRUCT(MKPOLS((W,EW)))
	WW = AA(LIST)(range(len(W)))
	VIEW(larModelNumbering(1,1,1)(W,[WW,EW,FW],submodel,0.6)) 
	

	## Computation of oriented ∂_3

	B_3 = Boundary3(W,EW,FW,boundaryFaces)
	print B_3.todense()
	boundary = B_3 * totalChain(range(len(CW)))
	b = boundary.tocoo()
	SB_2 = SBoundary2(EW,boundaryFaces)
	SB_2.todense()
	SB_2 = SB_2.tocsc()

	## Generation of oriented boundary faces
	
	polygons = []
	for face,sign in zip(b.row, b.data):
		print face,sign
		edges = list(SB_2[:,face].tocoo().row)
		signs = list(SB_2[:,face].tocoo().data)
		permutationMap = dict([EW[e] if sign*s>0 else REVERSE(EW[e]) for (e,s) in zip(edges,signs)])
		polygons += [orbit[:-1] for orbit in permutationOrbits(permutationMap)]
	# VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,polygons+EW))))	

	## Generation of oriented boundary triangles

	triangles = []
	for face in polygons:
		 triangles += AA(C(AL)(face[0]))(TRANS([face[1:-1],face[2:]]))
	P = W,triangles
	VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(P)))	

	## Computation of integral properties

	print "\nVolume(P) =", Volume(P)
	print "\nSurface(P) =", Surface(P)
