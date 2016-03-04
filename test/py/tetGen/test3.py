Extraction of 3-cells from a 2-skeleton embedded in 3D
========================================

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
    cycle = edgecyle + [edgecyle[0]]
    assert cycle[k] == 0
    if cycle[k+1] > 0: k_end = EV[k+1][0]
    else: k_end = EV[k+1][1]
    u,v = EV[k]
    if v==k_end: return True
    else: return False
    


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
"""
repeat until the set of non-worked faces is empty
"""
while nonWorkedFaces != set():
    """
    Take an elementary 2-chain F = [f]
    """
    seedFace = nonWorkedFaces.pop()
    nonWorkedFaces = nonWorkedFaces.difference({seedFace})
    faceChain = {seedFace}
    """
    compute the oriented 1-boundary  E := ∂_2(F) (loops of signed edges)
    """
    vect = csc_matrix((m,1),dtype='b')
    for face in faceChain:  vect[face] = 1
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
        for k,edge in enumerate(edgeCycle):
            print "\nedge=",edge
            """
            choose the "next" face g_i  on "ordered" coboundary of edge
            """
            if edge > 0:  faceLoop = REVERSE(EF_angle[edge])
            elif edge < 0:  faceLoop = EF_angle[-edge]
            elif edge == 0: 
                if zeroPos(EV,edgeCycle,k): 
                    faceLoop = REVERSE(EF_angle[edge])
                else:
                    faceLoop = EF_angle[-edge]
            faceLoop = faceLoop + [faceLoop[0]]
            print "faceLoop=",faceLoop
            pivot = faceChain.intersection(faceLoop).pop()
            print "pivot=",pivot
            pivotIndex = faceLoop.index(pivot)
            print "pivotIndex=",pivotIndex
            adjacentFace = faceLoop[pivotIndex+1]
            """
            assemble all the g_i with F in a new 2-chain F := F \cup [g_i]
            """
            faceChain = faceChain.union([adjacentFace])
            print "faceChain=",faceChain
            #VIEW(STRUCT(MKTRIANGLES((V,[FV[f] for f in faceChain],EV),color=True)))
            nonWorkedFaces = nonWorkedFaces.difference([adjacentFace])
            print "nonWorkedFaces=",nonWorkedFaces
            """
            """
        vect = csc_matrix((m,1),dtype='b')
        for face in faceChain:  vect[face] = 1
        edgeCycleCoords = boundaryOperator * vect
        edgeCycle = coords2chain(edgeCycleCoords)
        print "\nedgeCycle=",edgeCycle
    """
    put the signed F elements in a new column of ∂_3 (new row of coboundary_2)
    """

        
        