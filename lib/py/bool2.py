""" Module for Boolean computations between geometric objects """
from pyplasm import *
""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')
from inters import *
DEBUG = False

""" Coding utilities """
global count
""" Generation of a random 3D point """
def rpoint():
    return eval( vcode([ random.random(), random.random(), random.random() ]) )

""" Generation of random triangles """
def randomTriangles(numberOfTriangles=400,scaling=0.3):
    randomTriaArray = [rtriangle(scaling) for k in range(numberOfTriangles)]
    [xs,ys,zs] = TRANS(CAT(randomTriaArray))
    xmin, ymin, zmin = min(xs), min(ys), min(zs)
    v = array([-xmin,-ymin, -zmin])
    randomTriaArray = [[list(v1+v), list(v2+v), list(v3+v)] for v1,v2,v3 in randomTriaArray]
    return randomTriaArray

""" Generation of random 3D quadrilaterals """
def randomQuads(numberOfQuads=400,scaling=0.3):
    randomTriaArray = [rtriangle(scaling) for k in range(numberOfQuads)]
    [xs,ys,zs] = TRANS(CAT(randomTriaArray))
    xmin, ymin, zmin = min(xs), min(ys), min(zs)
    v = array([-xmin,-ymin, -zmin])
    randomQuadArray = [AA(list)([ v1+v, v2+v, v3+v, v+v2-v1+v3 ]) for v1,v2,v3 in randomTriaArray]
    return randomQuadArray

""" Generation of a single random triangle """
def rtriangle(scaling):
    v1,v2,v3 = array(rpoint()), array(rpoint()), array(rpoint())
    c = (v1+v2+v3)/3
    pos = rpoint()
    v1 = (v1-c)*scaling + pos
    v2 = (v2-c)*scaling + pos
    v3 = (v3-c)*scaling + pos
    return tuple(eval(vcode(v1))), tuple(eval(vcode(v2))), tuple(eval(vcode(v3)))

""" Containment boxes """
def containmentBoxes(randomPointArray,qualifier=0):
    if len(randomPointArray[0])==2:
        boxes = [eval(vcode([min(x1,x2), min(y1,y2), min(z1,z2), 
                             max(x1,x2), max(y1,y2), max(z1,z2)]))+[qualifier]
                for ((x1,y1,z1),(x2,y2,z2)) in randomPointArray]
    elif len(randomPointArray[0])==3:
        boxes = [eval(vcode([min(x1,x2,x3), min(y1,y2,y3), min(z1,z2,z3), 
                             max(x1,x2,x3), max(y1,y2,y3), max(z1,z2,z3)]))+[qualifier]
                for ((x1,y1,z1),(x2,y2,z2),(x3,y3,z3)) in randomPointArray]
    elif len(randomPointArray[0])==4:
        boxes = [eval(vcode([min(x1,x2,x3,x4), min(y1,y2,y3,y4), min(z1,z2,z3,z4), 
                             max(x1,x2,x3,x4), max(y1,y2,y3,y4), max(z1,z2,z3,z4)]))+[qualifier]
                for ((x1,y1,z1),(x2,y2,z2),(x3,y3,z3),(x4,y4,z4)) in randomPointArray]
    return boxes

""" Transformation of a 3D box into an hexahedron """    
def box2exa(box):
    x1,y1,z1,x2,y2,z2,type = box
    verts = [[x1,y1,z1], [x1,y1,z2], [x1,y2,z1], [x1,y2,z2], [x2,y1,z1], [x2,y1,z2], [x2,y2,z1], [x2,y2,z2]]
    cell = [range(1,len(verts)+1)]
    return [verts,cell,None],type

def lar2boxes(model,qualifier=0):
    V,CV = model
    boxes = []
    for k,cell in enumerate(CV):
        verts = [V[v] for v in cell]
        x1,y1,z1 = [min(coord) for coord in TRANS(verts)]
        x2,y2,z2 = [max(coord) for coord in TRANS(verts)]
        boxes += [eval(vcode([min(x1,x2),min(y1,y2),min(z1,z2),max(x1,x2),max(y1,y2),max(z1,z2)]))+[(qualifier,k)]]
    return boxes

""" Computation of the 1D centroid of a list of 3D boxes """    
def centroid(boxes,coord):
    delta,n = 0,len(boxes)
    ncoords = len(boxes[0])/2
    a = coord%ncoords
    b = a+ncoords
    for box in boxes:
        delta += (box[a] + box[b])/2
    return delta/n

""" Characteristic matrix transposition """
def invertRelation(CV):
    def myMax(List):
        if List==[]:  return -1
        else:  return max(List)
            
    columnNumber = max(AA(myMax)(CV))+1
    VC = [[] for k in range(columnNumber)]
    for k,cell in enumerate(CV):
        for v in cell: VC[v] += [k]
    return VC

""" Generation of a list of HPCs from a LAR model with non-convex faces """
def MKTRIANGLES(*model): 
    V,FV = model
    triangleSets = boundaryTriangulation(V,FV)
    return [ STRUCT([MKPOL([verts,[[1,2,3]],None]) for verts in triangledFace]) 
        for triangledFace in triangleSets ]

def MKSOLID(*model): 
    V,FV = model
    pivot = V[0]
    VF = invertRelation(FV) 
    faces = [face for face in FV if face not in VF[0]]
    triangleSets = boundaryTriangulation(V,faces)
    return XOR([ MKPOL([face+[pivot], [range(1,len(face)+2)],None])
        for face in CAT(triangleSets) ])


""" Split the boxes between the below,above subsets """
def splitOnThreshold(boxes,subset,coord):
    theBoxes = [boxes[k] for k in subset]
    threshold = centroid(theBoxes,coord)
    ncoords = len(boxes[0])/2
    a = coord%ncoords
    b = a+ncoords
    below,above = [],[]
    for k in subset:
        if boxes[k][a] <= threshold: below += [k]
    for k in subset:
        if boxes[k][b] >= threshold: above += [k]
    return below,above

""" Test if bucket OK or append to splitting stack """
def splitting(bucket,below,above, finalBuckets,splittingStack):
    if (len(below)<4 and len(above)<4) or len(set(bucket).difference(below))<7 \
        or len(set(bucket).difference(above))<7: 
        finalBuckets.append(below)
        finalBuckets.append(above)
    else: 
        splittingStack.append(below)
        splittingStack.append(above)


""" Remove subsets from bucket list """
def removeSubsets(buckets):
    n = len(buckets)
    A = zeros((n,n))
    for i,bucket in enumerate(buckets):
        for j,bucket1 in enumerate(buckets):
            if set(bucket).issubset(set(bucket1)):
                A[i,j] = 1
    B = AA(sum)(A.tolist())
    out = [bucket for i,bucket in enumerate(buckets) if B[i]==1]
    return out

def geomPartitionate(boxes,buckets):
    geomInters = [set() for h in range(len(boxes))]
    for bucket in buckets:
        for k in bucket:
            geomInters[k] = geomInters[k].union(bucket)
    for h,inters in enumerate(geomInters):
        geomInters[h] = geomInters[h].difference([h])
    return AA(list)(geomInters)

""" Iterate the splitting until \texttt{splittingStack} is empty """
def boxTest(boxes,h,k):
    B1,B2,B3,B4,B5,B6,_ = boxes[k]
    b1,b2,b3,b4,b5,b6,_ = boxes[h]
    return not (b4<B1 or B4<b1 or b5<B2 or B5<b2 or b6<B3 or B6<b3)

def boxBuckets(boxes):
    bucket = range(len(boxes))
    splittingStack = [bucket]
    finalBuckets = []
    while splittingStack != []:
        bucket = splittingStack.pop()
        below,above = splitOnThreshold(boxes,bucket,1)
        below1,above1 = splitOnThreshold(boxes,above,2)
        below2,above2 = splitOnThreshold(boxes,below,2) 
               
        below11,above11 = splitOnThreshold(boxes,above1,3)
        below21,above21 = splitOnThreshold(boxes,below1,3)        
        below12,above12 = splitOnThreshold(boxes,above2,3)
        below22,above22 = splitOnThreshold(boxes,below2,3)  
              
        splitting(above1,below11,above11, finalBuckets,splittingStack)
        splitting(below1,below21,above21, finalBuckets,splittingStack)
        splitting(above2,below12,above12, finalBuckets,splittingStack)
        splitting(below2,below22,above22, finalBuckets,splittingStack)
        
        finalBuckets = list(set(AA(tuple)(finalBuckets)))
    parts = geomPartitionate(boxes,finalBuckets)
    parts = [[h for h in part if boxTest(boxes,h,k)] for k,part in enumerate(parts)]
    return AA(sorted)(parts)

""" Computation of affine face transformations """
def COVECTOR(points):
    pointdim = len(points[0])
    plane = Planef.bestFittingPlane(pointdim,
                    [item for sublist in points for item in sublist])
    return [plane.get(I) for I in range(0,pointdim+1)]

def faceTransformations(facet):
    covector = COVECTOR(facet)
    translVector = facet[0]
    # translation 
    newFacet = [ VECTDIFF([v,translVector]) for v in facet ]
    # linear transformation: boundaryFacet -> standard (d-1)-simplex
    d = len(facet[0])
    m = mat( newFacet[1:d] + [covector[1:]] )
    if det(m)==0.0:
        for k in range(len(facet)-2):
            m = mat( newFacet[1+k+1:d+k+1] + [covector[1:]] )
            if det(m)!=0.0: break
    transformMat = m.T.I
    # transformation in the subspace x_d = 0
    out = (transformMat * (mat(newFacet).T)).T.tolist()
    return transformMat


""" Computation of topological relation """
def crossRelation(XV,YV):
    csrXV = csrCreate(XV)
    csrYV = csrCreate(YV)
    csrXY = matrixProduct(csrXV, csrYV.T)
    XY = [None for k in range(len(XV))]
    for k,face in enumerate(XV):
        data = csrXY[k].data
        col = csrXY[k].indices
        XY[k] = [col[h] for h,val in enumerate(data) if val==2] 
        # NOTE: val depends on the relation under consideration ...
    return XY

""" Submanifold mapping computation """
def submanifoldMapping(pivotFace):
    tx,ty,tz = pivotFace[0]
    transl = mat([[1,0,0,-tx],[0,1,0,-ty],[0,0,1,-tz],[0,0,0,1]])
    facet = [ VECTDIFF([v,pivotFace[0]]) for v in pivotFace ]
    m = faceTransformations(facet)
    mapping = mat([[m[0,0],m[0,1],m[0,2],0],[m[1,0],m[1,1],m[1,2],0],[m[2,0],
                    m[2,1],m[2,2],0],[0,0,0,1]])
    transform = mapping * transl
    return transform

""" Set of line segments partitioning a facet """
def intersection(V,FV,EV):
    def intersection0(k,bundledFaces):
        FE = crossRelation(FV,EV)
        pivotFace = [V[v] for v in FV[k]]
        transform = submanifoldMapping(pivotFace)  # submanifold transformation
        transformedCells,edges,faces = [],[],[]
        for face in bundledFaces:
            edge = set(FE[k]).intersection(FE[face])  # common edge index
            if edge == set():
                candidateEdges = FE[face]
                facet = []
                for e in candidateEdges:
                    cell = [V[v]+[1.0] for v in EV[e]]  # verts of incident face
                    transformedCell = (transform * (mat(cell).T)).T.tolist()  
                    # vertices in local frame
                    facet += [[point[:-1] for point in transformedCell]]
                faces += [facet]
            else:  # boundary edges of face k
                e, = edge
                vs = [V[v]+[1.0] for v in EV[e]]
                ws = (transform * (mat(vs).T)).T.tolist()
                edges += [[p[:-1] for p in ws]]
        return edges,faces,transform
    return intersection0    

""" Space partitioning via submanifold mapping """
def spacePartition(V,FV,EV, parts):
    transfFaces = []
    for k,bundledFaces in enumerate(parts):
        edges,faces,transform = intersection(V,FV,EV)(k,bundledFaces)
        for face in faces:
            line = []
            for edge in face:
                (x1,y1,z1),(x2,y2,z2) = edge
                if not verySmall(z2-z1):
                    x = (x2-x1)/(z2-z1) + x1
                    y = (y2-y1)/(z2-z1) + y1
                    p = [x,y,0]
                    line += [eval(vcode(p))]
            if line!=[]: edges += [line]
        v,fv,ev = larFromLines([[point[:-1] for point in edge] for edge in edges])    
        if len(fv)>1: fv = fv[:-1]
        lar = [w+[0.0] for w in v],fv,ev
        transfFaces += [Struct([ larApply(transform.I)(lar) ])]
    W,FW,EW = struct2lar(Struct(transfFaces))
    return W,FW,EW

from support import PolygonTessellator,vertex

def orientTriangle(pointTriple):
    v1 = array(pointTriple[1])-pointTriple[0]
    v2 = array(pointTriple[2])-pointTriple[0]
    if cross(v1,v2)[2] < 0: return REVERSE(pointTriple)
    else: return pointTriple

def boundaryTriangulation(W,FW):
    triangleSet = []
    for face in FW:
        pivotFace = [W[v] for v in face+(face[0],)]
        transform = submanifoldMapping(pivotFace)
        mappedVerts = (transform * (mat([p+[1.0] for p in pivotFace]).T)).T.tolist()
        facet = [point[:-2] for point in mappedVerts]
        pol = PolygonTessellator()
        vertices = [ vertex.Vertex( (x,y,0) ) for (x,y) in facet  ]
        verts = pol.tessellate(vertices)
        ps = [list(v.point) for v in verts]
        trias = [[ps[k],ps[k+1],ps[k+2],ps[k]] for k in range(0,len(ps),3)]
        mappedVerts = (transform.I * (mat([p+[1.0] for p in ps]).T)).T.tolist()
        points = [p[:-1] for p in mappedVerts]
        trias = [[points[k],points[k+1],points[k+2],points[k]] 
            for k in range(0,len(points),3) 
            if scipy.linalg.norm(cross(array(points[k+1])-points[k], 
                                       array(points[k+2])-points[k])) != 0 ]
        triangleSet += [AA(orientTriangle)(trias)]
    return triangleSet


def triangleIndices(triangleSet,W):
    vertDict,out = defaultdict(),[]
    for k,vertex in enumerate(W):  vertDict[vcode(vertex)] = k
    for h,faceSetOfTriangles in enumerate(triangleSet):
        out += [[[vertDict[vcode(p)] for p in triangle[:-1]] 
                    for triangle in faceSetOfTriangles]]
    return out


def edgesTriangles(EF, FW, TW, EW):
    ET = [None for k in range(len(EF))]
    for e,edgeFaces in enumerate(EF):
        ET[e] = []
        for f in edgeFaces:
            for t in TW[f]:
                if set(EW[e]).intersection(t)==set(EW[e]):
                    ET[e] += [t]
    return ET

""" Directional and orthogonal projection operators """
def dirProject (e):
    def dirProject0 (v):
        return SCALARVECTPROD([ INNERPROD([ UNITVECT(e), v ]), UNITVECT(e) ])
    return dirProject0

def orthoProject (e):
    def orthoProject0 (v):  
        return VECTDIFF([ v, dirProject(UNITVECT(e))(v) ])
    return orthoProject0

""" Circular ordering of faces around edges """

def planeProjection(normals):
    V = mat(normals)
    if all(V[:,0]==0): V = np.delete(V, 0, 1)
    elif all(V[:,1]==0): V = np.delete(V, 1, 1)
    elif all(V[:,2]==0): V = np.delete(V, 2, 1)
    return V

def faceSlopeOrdering(model):
    V,FV,EV = model
    triangleSet = boundaryTriangulation(V,FV)
    TV = triangleIndices(triangleSet,V)
    triangleVertices = CAT(TV)
    TE = crossRelation(triangleVertices,EV)
    ET,ET_angle = invertRelation(TE),[]
    for e,et in enumerate(ET):
        v1,v2 = EV[e]
        v1v2 = set([v1,v2])
        et_angle = []
        t0 = et[0]
        tverts = [v1,v2] + list(set(triangleVertices[t0]).difference(v1v2))
        e3 = UNITVECT(VECTDIFF([ V[tverts[1]], V[tverts[0]] ]))
        e1 = UNITVECT(VECTDIFF([ V[tverts[2]], V[tverts[0]] ]))
        e2 = cross(array(e1),e3).tolist()
        basis = mat([e1,e2,e3]).T
        transform = basis.I
        normals = []
        Tvs = []
        for triangle in et:
            verts = triangleVertices[triangle]
            vertSet = set(verts).difference(v1v2)
            tvs = [v1,v2] + list(vertSet)
            Tvs += [tvs]
            w1 = UNITVECT(VECTDIFF([ V[tvs[2]], V[tvs[0]] ]))
            w2 = (transform * mat([w1]).T).T
            w3 = cross(array([0,0,1]),w2).tolist()[0]
            normals += [w3]
        normals = mat(normals)
        for k,t in enumerate(et):
            angle = math.atan2(normals[k,1],normals[k,0])
            et_angle += [angle]
        pairs = sorted(zip(et_angle,et,Tvs))
        sortedTrias = [pair[1] for pair in pairs]
        triasVerts = [pair[2] for pair in pairs]
        #print "triasVerts =",triasVerts
        tetraVerts = triasVerts[0]+[triasVerts[1][2]]
        print det(mat([V[v]+[1] for v in tetraVerts]))
        #print "tetraVerts =",tetraVerts
        ET_angle += [sortedTrias]
    EF_angle = ET_to_EF_incidence(TV,FV, ET_angle)
    return EF_angle

""" Oriented cycle of vertices from a 1-cycle of unoriented edges """
def theNext(FE,EF_angle,EV,cb,previous_cb,previousOrientedEdges,cf):
    previous_cb = cb
    def theNext0(previous_edge,face):
        cbe = copy.copy(cb)
        edges = list(set(FE[face]).intersection(cbe)) #difference(cbe))
        if edges==[]: 
            edges = list(cbe)
            face = list(set(EF_angle[edges[0]]).intersection(cf))[0]
        if type(previousOrientedEdges[0])!=list:
            signs,next = cycles2permutation([previousOrientedEdges])
        else: signs,next = cycles2permutation(previousOrientedEdges)
        edge = edges[0]
        edgeOrientation = signs[edge]
        edgeFaces = EF_angle[edge]
        n = len(edgeFaces)
        if edgeOrientation == 1: 
            ind = (edgeFaces.index(face) + 1)%n
        elif edgeOrientation == -1:
            ind = (edgeFaces.index(face) - 1)%n
        nextFace = edgeFaces[ind]
        nextFaceBoundary = list(set(FE[nextFace]))
        orientedEdges = cyclesOrientation(previousOrientedEdges,nextFaceBoundary,EV)
        return orientedEdges,nextFace,edge
    return theNext0

""" Edge-triangles to Edge-faces incidence """
def ET_to_EF_incidence(TW,FW, ET_angle):
    tableFT = [None for k in range(len(FW))]
    t = 0
    for f,trias in enumerate(TW):
        tableFT[f] = range(t,t+len(trias))
        t += len(trias)
    tableTF = invertRelation(tableFT)
    EF_angle = [[tableTF[t][0] for t in triangles] for triangles in ET_angle]
    #assert( len(EF_angle) == 2*len(FW) )
    return EF_angle

""" Cells from $(d-1)$-dimensional LAR model """

def facesFromComponents(model):
    # initialization
    V,FV,EV = model
    EV = [EV[0]]+EV
    model = V,FV,EV
    EF_angle = faceSlopeOrdering(model)
    FE = crossRelation(FV,EV)
    # remove 0 indices from FE relation
    FE = [[f for f in face if f!=0] for face in FE]
    visitedCell = [[ None, None ] for k in range(len(FV)) ]
    print "\n>> 1: visitedCell =",[[k,row] for k,row in enumerate(visitedCell)]
    face = 0
    boundaryLoop = boundaryCycles(FE[face],EV)[0]
    firstEdge = boundaryLoop[0]
    cf,coe = getSolidCell(FE,face,visitedCell,boundaryLoop,EV,EF_angle,V,FV)
    for face,edge in zip(cf,coe):
        if visitedCell[face][0]==None: visitedCell[face][0] = edge
        else: visitedCell[face][1] = edge
    print "cf = ",cf
    print "coe = ",coe
    cv,ce = set(),set()
    cv = cv.union(CAT([FV[f] for f in cf]))
    ce = ce.union(CAT([FE[f] for f in cf]))
    CF,CV,CE,COE = [cf],[list(cv)],[list(ce)],[coe]
    
    # main loop
    while True:
        face, edge = startCell(visitedCell,FE,EV)
        if face == -1: break
        boundaryLoop = boundaryCycles(FE[face],EV)[0]
        if edge not in boundaryLoop:
            boundaryLoop = reverseOrientation(boundaryLoop)
        cf,coe = getSolidCell(FE,face,visitedCell,boundaryLoop,EV,EF_angle,V,FV)
        print "cf = ",cf
        print "coe = ",coe
        CF += [cf]
        COE += [coe]
        for face,edge in zip(cf,coe):
            if visitedCell[face][0]==None: visitedCell[face][0] = edge
            else: visitedCell[face][1] = edge
            
        cv,ce = set(),set()
        cv = cv.union(CAT([FV[f] for f in cf]))
        ce = ce.union(CAT([FE[f] for f in cf]))
        CV += [list(cv)]
        CE += [list(ce)]
    return V,CV,FV,EV,CF,CE,COE

""" Edge cycles associated to a closed chain of edges """
def boundaryCycles(edgeBoundary,EV):
    verts2edges = defaultdict(list)
    for e in edgeBoundary:
        verts2edges[EV[e][0]] += [e]
        verts2edges[EV[e][1]] += [e]
    cycles = []
    cbe = copy.copy(edgeBoundary)
    while cbe != []:
        e = cbe[0]
        v = EV[e][0]
        cycle = []
        while True:
            cycle += [(e,v)]
            e = list(set(verts2edges[v]).difference([e]))[0]
            cbe.remove(e)
            v = list(set(EV[e]).difference([v]))[0]
            if (e,v)==cycle[0]:
                break
        n = len(cycle)
        cycles += [[e if EV[e]==(cycle[(k-1)%n][1],cycle[k%n][1]) else -e 
            for k,(e,v) in enumerate(cycle)]]
    return cycles

""" Permutation of edges defined by edge cycles """
def cycles2permutation(cycles):
    next = []
    for cycle in cycles:
        next += zip(AA(ABS)(cycle),AA(ABS)(cycle[1:]+[cycle[0]]))
    next = dict(next)
    sign = dict([[ABS(edge),SIGN(edge)] for cycle in cycles for edge in cycle])
    return sign,next

""" Cycles orientation """
def cyclesOrient(pcycles,fcycle,EV):
    if set(AA(ABS)(pcycles)).difference(fcycle)==set(): return []
    ofcycle = boundaryCycles(fcycle,EV)[0] # oriented 
    if type(pcycles[0])==list: opcycle = CAT(pcycles)
    else: opcycle = pcycles
    int = set(opcycle).intersection(ofcycle)
    if int != set(): 
        ofcycle = reverseOrientation(ofcycle)
    outChain = [e for e in ofcycle if not (-e in opcycle)] 
    outChain += [e for e in opcycle if not (-e in ofcycle)] 
    return outChain

if __name__ == "__main__":
    pcycles = [[-19, 13, 22, 23]]
    fcycle = [30, 20, 18, 2, 26, 19]
    cyclesOrientation(pcycles,fcycle)

""" Start a new 3-cell """
def startCell(visitedCell,FE,EV):
    if len([term for cell in visitedCell for term in cell if term==None])==1: return -1,-1
    print "\n>> 3: visitedCell =",[[k,row] for k,row in enumerate(visitedCell)]
    for face in range(len(visitedCell)):
        if len([term for term in visitedCell[face] if term==None])==1:
            edge = visitedCell[face][0]
            print "\nface,edge =",face,edge
            break
        face,edge = -1,-1
    return face,edge

""" Face orientations storage """
def reverseOrientation(chain):
    return REVERSE([-cell for cell in chain])

def faceOrientation(boundaryLoop,face,FE,EV,cf):
    theBoundary = set(AA(ABS)(boundaryLoop))
    if theBoundary.intersection(FE[face])==set() and theBoundary.difference(FE[face])!=set(): ##BOH!!
        coboundaryFaces = [f for f in cf if set(FE[f]).intersection(theBoundary)!=set()]
        face = coboundaryFaces[0]            
    faceLoop = boundaryCycles(FE[face],EV)[0]
    commonEdges = set(faceLoop).intersection(boundaryLoop)
    if commonEdges == set() or commonEdges == {0}: 
        faceLoop = reverseOrientation(faceLoop)
        commonEdges = set(faceLoop).intersection(boundaryLoop)
    theEdge = list(commonEdges)[0]
    #if theEdge==0: theEdge = list(commonEdges)[1]
    return -theEdge,face

""" Check and store the orientation of faces """
def checkOrientation(previousOrientedEdges,orientedEdges,orientedFaceEdges,faceOrientations,face):
    list2 = CAT(orientedFaceEdges)
    if orientedEdges != []:
        list1 = CAT(orientedEdges)
    else: list1 = CAT(previousOrientedEdges)
    theList = set(list1).intersection(set(list2).union((lambda args:[-arg for arg in args])(list2)))
    if theList==set() or orientedEdges==[]:
        theList = set(CAT(orientedFaceEdges))
    edge = list(theList)[0]
    if theList.issubset(list1):  # equal signs
        if faceOrientations[face][0] == None:
            faceOrientations[face][0] = edge
        elif faceOrientations[face][1] == None:
            faceOrientations[face][1] = edge
        else: print "error: faceOrientations"
    elif not theList.issubset(list1): # different signs
        if faceOrientations[face][0] == None: 
            faceOrientations[face][0] = -edge
        elif faceOrientations[face][1] == None:
            faceOrientations[face][1] = -edge
        else: print "error: faceOrientations"
    else: print "error: checkOrientation"
    return faceOrientations

""" Get single solid cell """
def getSolidCell(FE,face,visitedCell,boundaryLoop,EV,EF_angle,V,FV):

    def orientFace(face,boundaryLoop): 
        for e in boundaryLoop:
            if ABS(e) in FE[face]: return -e
            
    coe = [orientFace(face,boundaryLoop)]
    cf = [face] 
    while boundaryLoop != []:
        edge,face = faceOrientation(boundaryLoop,face,FE,EV,cf)
        print "face,edge",face,edge
        if edge > 0: edgeFaces = EF_angle[edge]
        elif edge < 0: edgeFaces = REVERSE(EF_angle[-edge])
        e = ABS(edge)
        n = len(edgeFaces)
        ind = (edgeFaces.index(face)+1)%n
        nextFace = edgeFaces[ind]
        print "\nnextFace =",nextFace
        coe += [-orientFace(nextFace,boundaryLoop)]
        boundaryLoop = cyclesOrient(boundaryLoop,FE[nextFace],EV)
        print "boundaryLoop =",boundaryLoop
        cf += [nextFace] 
        face = nextFace
    print "\n>> 5: visitedCell =",[[k,row] for k,row in enumerate(visitedCell)]
    if DEBUG:
        VIEW(EXPLODE(1.2,1.2,1.2)( MKTRIANGLES(V,[FV[f] for f in cf]) ))
    return cf,coe

