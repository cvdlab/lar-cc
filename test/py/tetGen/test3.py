
"""
Extraction of 3-cells from a 2-skeleton embedded in 3D
========================================
"""

def coords2chain(chainCoords):
    chain = []
    chainCoordList = CAT(chainCoords.todense().tolist())
    for k,val in enumerate(chainCoordList):
        if val == 1:  chain += [k]
        elif val == -1:  chain += [-k]
    return chain

"""
compute the signed 2-boundary matrix, orienting every column (face) coherently with the orientation of its edge of greatest index (in order to remove the orientation indecision about edge 0 --- serve? Boh!)
"""
def larBoundary2(FV,EV):
    efOp = larFaces2Edges(FV,EV)
    FE = [efOp([k]) for k in range(len(FV))]
    data,row,col = [],[],[]
    for f in range(len(FE)):
        cycle_data = find_cycles(V,[EV[e] for e in FE[f]])
        coefficients = [cycle_data['ev_mapping'][k]['direction'] 
              for k in range(len(cycle_data['ev_mapping']))]
        ecycle = FE[f]
        data += coefficients
        row += ecycle
        col += [f]*len(ecycle)
        #print f,len(data),len(row),len(col)
    signedBoundary2 = coo_matrix((data,(row,col)), shape=(len(EV),len(FV)),dtype='b')
    return csr_matrix(signedBoundary2)

def zeroPos(EV,edgeCycle,k):
    cycle = edgeCycle + [edgeCycle[0]]
    assert cycle[k] == 0
    if cycle[k+1] > 0: k_end = EV[k+1][0]
    else: k_end = EV[k+1][1]
    u,v = EV[k]
    if v==k_end: return True
    else: return False

def zeroPos(EV,edgeCycle,k):
    cycle = edgeCycle + [edgeCycle[0]]
    assert cycle[k] == 0
    u,v = EV[k+1]
    if v>u and cycle[k+1] > 0: return True
    elif v>u and cycle[k+1] < 0: return False
    else: print "error: orientation of edge 0"
    

"""
choose the "next" face g_i  on "ordered" coboundary of edge
"""
def adjFace(boundaryOperator,EV,EF_angle,faceChain,edgeCycle):
    def adjFace0(k,edge):
        if edge > 0:  edgeLoop = REVERSE(EF_angle[edge])
        elif edge < 0:  edgeLoop = EF_angle[-edge]
        elif edge == 0: 
            if zeroPos(EV,edgeCycle,k): 
                edgeLoop = REVERSE(EF_angle[edge])
            else:
                edgeLoop = EF_angle[-edge]
        edgeLoop = edgeLoop + [edgeLoop[0]]  # all positive indices
        print "edgeLoop =",edgeLoop
        pivotFace = set([ABS(f) for f in faceChain]).intersection(edgeLoop).pop()
        if pivotFace in edgeLoop:
            pivotIndex = edgeLoop.index(pivotFace)
        else:
            pivotIndex = edgeLoop.index(-pivotFace)
        adjacentFace = edgeLoop[pivotIndex+1]
        theSign = boundaryOperator[ABS(edge),adjacentFace]
        print "adjacentFace,edge,pivotFace,theSign =",adjacentFace,edge,pivotFace,theSign
        return adjacentFace * -(theSign*SIGN(edge))
    return adjFace0

def chooseSeedFace(FV,nonWorkedFaces,faceChain):
    if len(FV) == len(nonWorkedFaces): return nonWorkedFaces.pop()
    else: pass

def larBoundary3((V,FV,EV)):
	"""
	sort on angles the co-boundaries of 1-cells (loops of signed faces)
	"""
	model = V,FV,EV
	efOp = larFaces2Edges(FV,EV)
	FE = [efOp([k]) for k in range(len(FV))]
	EF_angle, _,_,_ = faceSlopeOrdering(model,FE)
	"""
	initialize the 2-coboundary array of arrays, with boundary 2-faces of 3-cells by row
	"""
	m = len(FV)
	nonWorkedFaces = set(range(m))
	coboundary_2 = []
	boundaryOperator = larBoundary2(FV,EV)
	cellNumber=0
	row,col,data = [],[],[]
	faceChain,CF = {},[]
	"""
	repeat until the set of non-worked faces is empty
	"""
	while nonWorkedFaces != set():
		"""
		Take an elementary 2-chain F = [f]
		"""
		seedFace = nonWorkedFaces.pop()
		#seedFace = chooseSeedFace(FV,nonWorkedFaces)
		#nonWorkedFaces = nonWorkedFaces.difference({seedFace})
		faceChain = {seedFace}
		"""
		compute the oriented 1-boundary  E := ∂_2(F) (loops of signed edges)
		"""
		vect = csc_matrix((m,1),dtype='b')
		for face in faceChain:  vect[face] = SIGN(face)
		edgeCycleCoords = boundaryOperator * vect
		edgeCycle = coords2chain(edgeCycleCoords)
		print "\nedgeCycle=",edgeCycle
		"""
		repeat until E[F] = ø
		"""
		while edgeCycle != []:
			"""
			for each edge on E 
			"""
			look4face = adjFace(boundaryOperator,EV,EF_angle,faceChain,edgeCycle)
			for k,edge in enumerate(edgeCycle):
				print "\nedge=",edge
				adjacentFace = look4face(k,edge)
				"""
				assemble all the g_i with F in a new 2-chain F := F \cup [g_i]
				"""
				faceChain = faceChain.union([adjacentFace])
				print "faceChain=",faceChain
				#VIEW(STRUCT(MKTRIANGLES((V,[FV[f] for f in faceChain],EV),color=True)))
				#nonWorkedFaces = nonWorkedFaces.difference([ABS(adjacentFace)])
				#print "nonWorkedFaces=",nonWorkedFaces
				"""
				"""
			vect = csc_matrix((m,1),dtype='b')
			for face in faceChain:  vect[ABS(face)] = SIGN(face)
			edgeCycleCoords = boundaryOperator * vect
			edgeCycle = coords2chain(edgeCycleCoords)
			print "\nedgeCycle=",edgeCycle
		"""
		put the signed F elements in a new column of ∂_3 (new row of coboundary_2)
		"""
		row += [ABS(face) for face in faceChain]
		col += [cellNumber for face in faceChain]
		data += [SIGN(face) for face in faceChain]
		cellNumber += 1
		CF += [list(faceChain)] 
		nonWorkedFaces = nonWorkedFaces.difference(AA(ABS)(faceChain))
		print "nonWorkedFaces=",nonWorkedFaces

	outMatrix = coo_matrix((data, (row,col)), shape=(m,cellNumber),dtype='b')
	CV = [list(set(CAT([FV[f] for f in cell]))) for cell in CF[:3]]
	return csr_matrix(outMatrix),CV


csrboundary3,CV = larBoundary3((V,FV,EV))

"""
V,BF,BE = larBoundary3(V,FV,EV,VV)([1]*len(FV))
V,BF,BE = larBoundary3(V,FV,EV,VV)([0]*3 +[1] +[0]*8 +[1])


VIEW(STRUCT(MKTRIANGLES((V,BF,BE),color=True))) 
VIEW(EXPLODE(1.2,1.2,1.2)(MKTRIANGLES((V,BF,BE),color=True))) 
"""