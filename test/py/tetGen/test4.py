

\subsection{Extraction of 3-cells from a 2-skeleton embedded in 3D}
Need to characterize the two vertices of zero edge as "last" of previous edge in the chain,
and "first" of following edge in the chain. Therefore the zero edge is oriented as from "last" to "first" and consequently, is "positive" iff  "first" > "last".

It is sufficient to characterize the "first" vertex of the edge following the zero edge. This one is the "last" of zero edge. the other vertex of zero edge is its "first". Zero edge is "positive" iff last(zero) > first(zero).

\paragraph{aaaaaaa}
%-------------------------------------------------------------------------------
@D aaaaaaa
@{""" aaaaaaa """
def coords2chain(chainCoords):
    coo = coo_matrix(chainCoords)
    return [(e,val) for e,val in zip(coo.row,coo.data)]
@}
%-------------------------------------------------------------------------------


\paragraph{aaaaaaa}
choose the "next" face g_i  on "ordered" coboundary of edge

%-------------------------------------------------------------------------------
@D aaaaaaa
@{""" aaaaaaa """
def adjFace(boundaryOperator,EV,EF_angle,faceChainOrientation):
    def adjFace0(edge,orientation):
        if orientation > 0:  edgeLoop = REVERSE(EF_angle[edge])
        elif orientation < 0:  edgeLoop = EF_angle[edge]
        edgeLoop = edgeLoop + [edgeLoop[0]]  # all positive indices
        print "edgeLoop =",edgeLoop
        
        print "faceChainOrientation =",faceChainOrientation
        pivotFace = set([f for f,_ in faceChainOrientation]).intersection(edgeLoop).pop()
        if pivotFace in edgeLoop:
            pivotIndex = edgeLoop.index(pivotFace)
        else:
            pivotIndex = edgeLoop.index(-pivotFace)
        adjacentFace = edgeLoop[pivotIndex+1]
        
        theSign = boundaryOperator[edge,adjacentFace]
        print "adjacentFace,edge,pivotFace,theSign =",adjacentFace,edge,pivotFace,theSign
        return adjacentFace, -(theSign*orientation)
        #return adjacentFace, theSign
    return adjFace0
@}
%-------------------------------------------------------------------------------


\paragraph{aaaaaaa}
%-------------------------------------------------------------------------------
@D aaaaaaa
@{""" aaaaaaa """
def chooseStartFace(FV,faceCounter):
    print "\n>>>>>>> ECCOMI"
    for f in range(len(FV)):
        if faceCounter[f,0]==1 and faceCounter[f,1]==0: return (f,-1)
        elif faceCounter[f,0]==0 and faceCounter[f,1]==1: return (f,1)
    if sum(array(faceCounter))==2*len(FV): return (-1,999)
    else: return (0,1)
@}
%-------------------------------------------------------------------------------


\paragraph{aaaaaaa}
%-------------------------------------------------------------------------------
@D aaaaaaa
@{""" aaaaaaa """
def signedBasis(boundaryOperator):
	facesByEdges = csc_matrix(boundaryOperator)
	m,n = facesByEdges.shape
	edges,signs = [],[]
	for i in range(n):
		edges += [facesByEdges.indices[facesByEdges.indptr[i]:facesByEdges.indptr[i+1]].tolist()]
		signs += [facesByEdges.data[facesByEdges.indptr[i]:facesByEdges.indptr[i+1]].tolist()]
	return zip(edges,signs)
@}
%-------------------------------------------------------------------------------



\paragraph{aaaaaaa}
%-------------------------------------------------------------------------------
@D aaaaaaa
@{""" aaaaaaa """
def larSignedBoundary3((V,FV,EV)):
	"""
	sort on angles the co-boundaries of 1-cells (loops of signed faces)
	"""
	model = V,FV,EV
	faceCounter = zeros((len(FV),2),dtype='b')
	CF = []
	efOp = larFaces2Edges(V,FV,EV)
	FE = [efOp([k]) for k in range(len(FV))]
	EF_angle, _,_,_ = faceSlopeOrdering(model,FE)
	"""
	initialize the 2-coboundary array of arrays, with boundary 2-faces of 3-cells by row
	"""
	m = len(FV)
	nonWorkedFaces = set(range(m))
	coboundary_2 = []
	boundaryOperator = larSignedBoundary2(V,FV,EV)
	FEbasis = signedBasis(boundaryOperator)
	cellNumber=0
	row,col,data = [],[],[]
	"""
	repeat until the set of non-worked faces is empty
	"""
	while True:
		print "\n>>>>>>> CIAO"
		"""
		Take an elementary 2-chain F = [f]
		"""
		#startFace = nonWorkedFaces.pop()
		startFace,orientation = chooseStartFace(FV,faceCounter)
		print "\n\n*** startFace =",startFace
		if startFace == -1: break
		nonWorkedFaces = nonWorkedFaces.difference({startFace})
		faceChainOrientation = {(startFace,orientation)}
		"""
		compute the oriented 1-boundary  E := ∂_2(F) (loops of signed edges)
		"""
		vect = csc_matrix((m,1),dtype='b')
		for face,orientation in faceChainOrientation:  
		    vect[face] = orientation
		edgeCycleCoords = boundaryOperator * vect
		edgeCycle = coords2chain(edgeCycleCoords)
		print "\nedgeCycle esterno=",edgeCycle
		"""
		repeat until E[F] = ø
		"""
		while edgeCycle != []:
			"""
			for each edge on E 
			"""
			look4face = adjFace(boundaryOperator,EV,EF_angle,faceChainOrientation)
			for edge,orientation in edgeCycle:
				print "\nedge,orientation=",edge,orientation
				adjacentFace,orientation = look4face(edge,orientation)
				"""
				assemble all the g_i with F in a new 2-chain F := F \cup [g_i]
				"""
				faceChainOrientation = faceChainOrientation.union(
				    [(adjacentFace,orientation)])
				print "faceChainOrientation=",faceChainOrientation
				#VIEW(STRUCT(MKTRIANGLES((V,[FV[f] for f in faceChain],EV),color=True)))
				nonWorkedFaces = nonWorkedFaces.difference([adjacentFace])
				print "nonWorkedFaces=",nonWorkedFaces
				"""
				"""
			vect = csc_matrix((m,1),dtype='b')
			for face,orientation in faceChainOrientation:  
			    vect[face] = orientation
			edgeCycleCoords = boundaryOperator * vect
			print "\nedgeCycle interno pre=",edgeCycle
			edgeCycle = coords2chain(edgeCycleCoords)
			print "\nedgeCycle interno post=",edgeCycle
		"""
		put the signed F elements in a new column of ∂_3 (new row of coboundary_2)
		"""
		row += [face for face,_ in faceChainOrientation]
		col += [cellNumber for face,orientation in faceChainOrientation]
		data += [orientation for _,orientation in faceChainOrientation]
		cellNumber += 1
		
		for face,orientation in faceChainOrientation:
		    if orientation == 1: faceCounter[face,0]+=1
		    elif orientation == -1: faceCounter[face,1]+=1
		print "faceChainOrientation =",faceChainOrientation	
		print "faceCounter =",faceCounter	
		    
		CF += [[face for face,_ in faceChainOrientation]]
		print "CF =",CF

	outMatrix = coo_matrix((data, (row,col)), shape=(m,cellNumber),dtype='b')
	return csr_matrix(outMatrix),CF,faceCounter
@}
%-------------------------------------------------------------------------------


\paragraph{aaaaaaa}
%-------------------------------------------------------------------------------
@D aaaaaaa
@{""" aaaaaaa """
if __name__=="__main__":

	V,[VV,EV,FV,CV] = larCuboids([2,1,1],True)
	cubeGrid = Struct([(V,FV,EV)],"cubeGrid")
	cubeGrids = Struct(2*[cubeGrid,t(.5,.5,.5),r(0,0,PI/6)])

	V,FV,EV = struct2Marshal(cubeGrids)
	csrmat,CF,faceCounter = larSignedBoundary3((V,FV,EV))
	csrmat.todense()
@}
%-------------------------------------------------------------------------------
