""" Module for Boolean computations between geometric objects """
from larlib import *
from copy import copy
DEBUG = False

""" Coding utilities """
global count
""" Generation of a random 3D point """
def rpoint3d():
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
    v1,v2,v3 = array(rpoint3d()), array(rpoint3d()), array(rpoint3d())
    c = (v1+v2+v3)/3
    pos = rpoint3d()
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

""" Generation of a list of HPCs from a LAR model with non-convex faces """
def MKTRIANGLES(model): 
    V,FV,EV = model
    if len(V[0]) == 2: V=[v+[0] for v in V]
    FE = crossRelation(V,FV,EV)
    triangleSets = boundaryTriangulation(V,FV,EV,FE)
    return [ STRUCT([MKPOL([verts,[[1,2,3]],None]) for verts in triangledFace]) 
        for triangledFace in triangleSets ]

def MKSOLID(*model): 
    V,FV,EV = model
    FE = crossRelation(V,FV,EV)
    pivot = V[0]
    VF = invertRelation(FV) 
    faces = [face for face in FV if face not in VF[0]]
    triangleSets = boundaryTriangulation(V,faces,EV,FE)
    return XOR([ MKPOL([face+[pivot], [range(1,len(face)+2)],None])
        for face in CAT(triangleSets) ])



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

""" Iterate the splitting until \texttt{splittingStack} is empty """
def boxTest(boxes,h,k):
    if len(boxes[0])==4:
        B1,B2,B3,B4 = boxes[k]
        b1,b2,b3,b4 = boxes[h]
        return not (b3<B1 or B3<b1 or b4<B2 or B4<b2)
    else:
        B1,B2,B3,B4,B5,B6,_ = boxes[k]
        b1,b2,b3,b4,b5,b6,_ = boxes[h]
        return not (b4<B1 or B4<b1 or b5<B2 or B5<b2 or b6<B3 or B6<b3)

def boxBuckets3d(boxes):
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
from numpy import array
from scipy.linalg import lu
from scipy.linalg.basic import det

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
    if abs(det(m))<0.0001:
        for k in range(len(facet)-2):
            m = mat( newFacet[1+k+1:d+k+1] + [covector[1:]] )
            if abs(det(m))>0.0001: break
    transformMat = m.T.I
    # transformation in the subspace x_d = 0
    out = (transformMat * (mat(newFacet).T)).T.tolist()
    return transformMat


""" Computation of topological relation """
def crossRelation(V,XV,YV):
    csrXV = csrCreate(XV,lenV=len(V))
    csrYV = csrCreate(YV,lenV=len(V))
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

""" Helper functions for spacePartition """
def submodel(V,FV,EV):
    FE = crossRelation(V,FV,EV)
    def submodel0(f,F):
        fE = list(set(FE[f] + CAT([FE[g] for g in F])))
        fF = [f]+F
        return fF,fE
    return submodel0

def meetZero( sW, (w1,w2) ):
    testValue = sW[w1][2] * sW[w2][2]
    if testValue > 10**-4: 
        return False
    else: return True

def segmentIntersection(p1,p2):
    (x1,y1,z1),(x2,y2,z2) = p1,p2
    if abs(z1-z2) != 0.0:
        alpha = z1/(z1-z2)
        x = x1+alpha*(x2-x1)
        y = y1+alpha*(y2-y1)
        return x,y,0.0
    else: return None

""" Sorting of points along a line """
def computeSegments(line):
    p0 = line[0]
    p1 = line[1]
    p1_p0 = VECTDIFF([p1,p0])
    h = sorted([(X_k,k) for k,X_k in enumerate(p1_p0) if X_k!=0.0], reverse=True)[0][1]
    params = [(p[h]-p0[h])/p1_p0[h] for p in line]
    sortedPoints = TRANS(sorted(zip(params,line)))[1]
    segments = TRANS([sortedPoints[:-1],sortedPoints[1:]])
    return segments

""" Remove intersection points external to a submanifold face """
def removeExternals(M,V,EV,fe, z,fz,ez):
    inputFace = Struct([(V,[EV[e] for e in fe])])
    v,ev = larApply(M)(struct2lar(inputFace))
    pol = [eval(vcode(w[:-1])) for w in v],ev
    out = []
    classify = pointInPolygonClassification(pol)
    for k,point in enumerate(z):
        if classify(point)=="p_out":  out += [k]
    fz = [f for f in fz if not any([v in out for v in f])]
    ez = [e for e in ez if not any([v in out for v in e])]
    return z,fz,ez

""" Space partitioning via submanifold mapping """
def spacePartition(V,FV,EV, parts):
    FE = crossRelation(V,FV,EV)
    submodel0 = submodel(V,FV,EV)
    out = []
    """ input: face index f; candidate incident faces F[f]; """
    for f,F in enumerate(parts):
        """ Selection of the LAR submodel S(f) := (V,FV,EV)(f) restricted to [f]+F[f] """    
        fF,fE = submodel0(f,F)
        subModel = Struct([(V,[FV[g] for g in fF],[EV[h] for h in fE])])
        sV,sFV,sEV = struct2lar(subModel)
        """ Computation of submanifold map M moving f to z=0 """
        pivotFace = [V[v] for v in FV[f]]
        M = submanifoldMapping(pivotFace)  # submanifold transformation
        """ Transformation of S(f) by M, giving S = (sW,sEW) := M(S(f)) """
        sW,sFW,sEW = larApply(M)((sV,sFV,sEV))
        """ filtering of EW edges traversing z=0, giving EZ edges and incident faces FZEZ """
        sFE = crossRelation(V,sFW,sEW)    
        edges = list(set([ e for k,face in enumerate(sFW)  for e in sFE[k] 
                    if meetZero(sW, sEW[e]) ]))
        edgesPerFace = [ [e for e in sFE[k] if meetZero(sW, sEW[e])] 
                    for k,face in enumerate(sFW) ]
        edges = list(set(CAT(edgesPerFace)))
        WW = AA(LIST)(range(len(sW)))
        wires = [sEW[e] for e in edges]
        """ for each face in FZEZ, computation of the aligned set of points p(z=0) """
        points = collections.OrderedDict()
        lines = [[sW[w1],sW[w2]] for w1,w2 in wires]
        for k,(p,q) in enumerate(lines): 
            point = segmentIntersection(p,q)
            if point != None: points[edges[k]] = point
        pointsPerFace = [set(face).intersection(points.keys()) for face in edgesPerFace]
        lines = [[points[e][:2] for e in face] for face in pointsPerFace]
        lines = [line for line in lines if line!=[]]
        vpoints = [[(vcode(point),k) for k,point in enumerate(line)] for line in lines]
        lines = [AA(eval)(dict(line).keys()) for line in vpoints]
        """ sorting of every aligned set FX, where X is the parameter along the intersection line """
        theLines = []
        for line in lines:
            if len(line)>2: 
                print '\a',"line =",line
                segments = computeSegments(line)
                theLines += segments
            else: theLines += [line]
        ### Check that every set FX has even cardinality
        """ Construction of the planar set FX,EX of faces and lines """
        z,fz,ez = larFromLines(theLines)
        """ Remove external vertices """
        z,fz,ez = removeExternals(M,V,EV,FE[f], z,fz,ez)
        w,fw,ew = larApply(M.I)(([v+[0.0] for v in z],fz,ez))
        out += [Struct([(w,fw,ew)])]
    return struct2lar(Struct(out))


""" Half-line crossing test """
def crossingTest(new,old,count,status):
    if status == 0:
        status = new
        count += 0.5
    else:
        if status == old: count += 0.5
        else: count -= 0.5
        status = 0

""" Tile codes computation """
def setTile(box):
    tiles = [[9,1,5],[8,0,4],[10,2,6]]
    b1,b2,b3,b4 = box
    def tileCode(point):
        x,y = point
        code = 0
        if y>b1: code=code|1
        if y<b2: code=code|2
        if x>b3: code=code|4
        if x<b4: code=code|8
        return code 
    return tileCode

""" Point in polygon classification """
def pointInPolygonClassification(pol):

    V,EV = pol
    # edge orientation
    FV = [sorted(set(CAT(EV)))]
    orientedCycles = boundaryPolylines(Struct([(V,FV,EV)]))
    EV = []
    for cycle in orientedCycles:
        EV += zip(cycle[:-1],cycle[1:])

    def pointInPolygonClassification0(p):
        x,y = p
        xmin,xmax,ymin,ymax = x,x,y,y
        tilecode = setTile([ymax,ymin,xmax,xmin])
        count,status = 0,0
    
        for k,edge in enumerate(EV):
            p1,p2 = edge[0],edge[1]
            (x1,y1),(x2,y2) = p1,p2
            c1,c2 = tilecode(p1),tilecode(p2)
            c_edge, c_un, c_int = c1^c2, c1|c2, c1&c2
            
            if c_edge == 0 and c_un == 0: return "p_on"
            elif c_edge == 12 and c_un == c_edge: return "p_on"
            elif c_edge == 3:
                if c_int == 0: return "p_on"
                elif c_int == 4: count += 1
            elif c_edge == 15:
                x_int = ((y-y2)*(x1-x2)/(y1-y2))+x2 
                if x_int > x: count += 1
                elif x_int == x: return "p_on"
            elif c_edge == 13 and ((c1==4) or (c2==4)):
                    crossingTest(1,2,status,count)
            elif c_edge == 14 and (c1==4) or (c2==4):
                    crossingTest(2,1,status,count)
            elif c_edge == 7: count += 1
            elif c_edge == 11: count = count
            elif c_edge == 1:
                if c_int == 0: return "p_on"
                elif c_int == 4: crossingTest(1,2,status,count)
            elif c_edge == 2:
                if c_int == 0: return "p_on"
                elif c_int == 4: crossingTest(2,1,status,count)
            elif c_edge == 4 and c_un == c_edge: return "p_on"
            elif c_edge == 8 and c_un == c_edge: return "p_on"
            elif c_edge == 5:
                if (c1==0) or (c2==0): return "p_on"
                else: crossingTest(1,2,status,count)
            elif c_edge == 6:
                if (c1==0) or (c2==0): return "p_on"
                else: crossingTest(2,1,status,count)
            elif c_edge == 9 and ((c1==0) or (c2==0)): return "p_on"
            elif c_edge == 10 and ((c1==0) or (c2==0)): return "p_on"
        if ((round(count)%2)==1): return "p_in"
        else: return "p_out"
    return pointInPolygonClassification0


""" From structures to boundary polylines """
def boundaryPolylines(struct):
    V,boundaryEdges = structBoundaryModel(struct)
    polylines = boundaryModel2polylines((V,boundaryEdges))
    return polylines

""" From LAR oriented boundary model to polylines """
def boundaryModel2polylines(model):
    V,EV = model
    polylines = []
    succDict = dict(EV)
    visited = [False for k in range(len(V))]
    nonVisited = [k for k in succDict.keys() if not visited[k]]
    while nonVisited != []:
        first = nonVisited[0]; v = first; polyline = []
        while visited[v] == False:
            visited[v] = True; 
            polyline += V[v], 
            v = succDict[v]
        polyline += [V[first]]
        polylines += [polyline]
        nonVisited = [k for k in succDict.keys() if not visited[k]]
    return polylines

""" 3D boundary triangulation of the space partition """

def orientTriangle(pointTriple):
    v1 = array(pointTriple[1])-pointTriple[0]
    v2 = array(pointTriple[2])-pointTriple[0]
    if cross(v1,v2)[2] < 0: return REVERSE(pointTriple)
    else: return pointTriple
    
from copy import copy

def boundaryTriangulation(V,FV,EV,FE):
    triangleSet = []  
    
    def mapVerts(inverseMap):
        def mapVerts0(mappedVerts):
            return (inverseMap * (mat(mappedVerts).T)).T.tolist()
        return mapVerts0
        
    for f,face in enumerate(FV):
        triangledFace = []
        EW = [EV[e] for e in FE[f]]
        pivotFace = [V[v] for v in face]
        vertdict = dict([(w,v) for v,w in enumerate(face)])
        EW = [[vertdict[w] for w in edge] for edge in EW]
        transform = submanifoldMapping(pivotFace)
        mappedVerts = (transform * (mat([p+[1.0] for p in pivotFace]).T)).T.tolist()
        verts2D = [point[:-2] for point in mappedVerts] 
              
        # reconstruction of boundary polyline for LAR face
        model = (verts2D,[range(len(verts2D))],EW)
        struct = Struct([model])
        points = boundaryModel2polylines(structBoundaryModel(struct))[0]
        
        # CDT triangulation with poly2tri
        polyline = [Point(p[0],p[1]) for p in points[:-1]]  
        cdt = CDT(REVERSE(polyline))
        triangles = cdt.triangulate()
        trias = [ [[t.c.x,t.c.y,0,1],[t.b.x,t.b.y,0,1],[t.a.x,t.a.y,0,1]] for t in triangles ]
        inverseMap = transform.I
        trias = AA(mapVerts(inverseMap))(trias)
        
        triangledFace += [[v[:-1] for v in triangle] for triangle in trias]
        triangleSet += [triangledFace]
    return triangleSet

def triangleIndices(triangleSet,W):
    vertDict,out = defaultdict(),[]
    for k,vertex in enumerate(W):  vertDict[vcode(vertex,PRECISION=3)] = k
    for h,faceSetOfTriangles in enumerate(triangleSet):
        trias = [[vertDict[vcode(p,PRECISION=3)] for p in triangle]
                    for triangle in faceSetOfTriangles]
        out += [trias]
    assert len(W)==max(CAT(CAT(out)))+1
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



""" Circular ordering of faces around edges """

def planeProjection(normals):
    V = mat(normals)
    if all(V[:,0]==0): V = np.delete(V, 0, 1)
    elif all(V[:,1]==0): V = np.delete(V, 1, 1)
    elif all(V[:,2]==0): V = np.delete(V, 2, 1)
    return V

def faceSlopeOrdering(model,FE,Z):
    V,FV,EV = model
    triangleSet = boundaryTriangulation(V,FV,EV,FE)
    TV = triangleIndices(triangleSet,Z)
    triangleVertices = CAT(TV)
    TE = crossRelation(V,triangleVertices,EV)
    ET,ET_angle = invertRelation(TE),[]
    #import pdb; pdb.set_trace()
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
        ET_angle += [sortedTrias]
    EF_angle = ET_to_EF_incidence(TV,FV, ET_angle)
    return EF_angle


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

def facesFromComponents(model,FE,EF_angle):
    # initialization
    V,FV,EV = model
    visitedCell = [[ None, None ] for k in range(len(FV)) ]
    face = 0
    boundaryLoop = boundaryCycles(FE[face],EV)[0]
    firstEdge = boundaryLoop[0]
    #import pdb; pdb.set_trace()
    cf,coe = getSolidCell(FE,face,visitedCell,boundaryLoop,EV,EF_angle,V,FV)
    for face,edge in zip(cf,coe):
        if visitedCell[face][0]==None: visitedCell[face][0] = edge
        else: visitedCell[face][1] = edge
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
    #cyclesOrientation(pcycles,fcycle)

""" Start a new 3-cell """
def startCell(visitedCell,FE,EV):
    if len([term for cell in visitedCell for term in cell if term==None])==1: return -1,-1
    for face in range(len(visitedCell)):
        if len([term for term in visitedCell[face] if term==None])==1:
            edge = visitedCell[face][0]
            break
        else: pass  #TODO: implement search for isolated shells
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


""" Get single solid cell """
def getSolidCell(FE,face,visitedCell,boundaryLoop,EV,EF_angle,V,FV):

    def orientFace(face,boundaryLoop): 
        for e in boundaryLoop:
            if ABS(e) in FE[face]: return -e

    coe = [orientFace(face,boundaryLoop)]
    cf = [face] 
    while boundaryLoop != []:
        edge,face = faceOrientation(boundaryLoop,face,FE,EV,cf)
        if edge > 0: edgeFaces = EF_angle[edge]
        elif edge < 0: edgeFaces = REVERSE(EF_angle[-edge])
        e = ABS(edge)
        n = len(edgeFaces)
        ind = (edgeFaces.index(face)+1)%n
        nextFace = edgeFaces[ind]
        coe += [-orientFace(nextFace,boundaryLoop)]
        boundaryLoop = cyclesOrient(boundaryLoop,FE[nextFace],EV)
        cf += [nextFace] 
        face = nextFace
    if DEBUG: pass
        #VIEW(EXPLODE(1.2,1.2,1.2)( MKTRIANGLES((V,[FV[f] for f in cf])) )) #add EV!
    return cf,coe

""" Main procedure of arrangement partitioning """

""" Double check the faces boundaries made of edges """
def doubleCheckFaceBoundaries(FE,V,FV,EV):
    FEout = []
    for f,face in enumerate(FE):
        n = len(FV[f])
        if len(FE[f]) > n:
            # contains both edges coded 0 and 1 ... (how to solve?)
            FEout += [list(set(face).difference([0]))]
        else:
            FEout += [face]
    return FEout


def thePartition(W,FW,EW):
    quadArray = [[W[v] for v in face] for face in FW]
    parts = boxBuckets3d(containmentBoxes(quadArray))
    Z,FZ,EZ = spacePartition(W,FW,EW, parts)
    Z,FZ,EZ = larSimplify((Z,FZ,EZ),radius=0.001)
    EZ = [EZ[0]]+EZ
    model = Z,FZ,EZ

    ZZ = AA(LIST)(range(len(Z)))
    """
    for k,face in enumerate(FZ):
        submodel = STRUCT(MKPOLS((Z,[face]+EZ)))
    """
    submodel = STRUCT(MKPOLS((Z,EZ)))
    VIEW(larModelNumbering(1,1,1)(Z,[ZZ,EZ,FZ],submodel,0.1)) 

    FE = crossRelation(Z,FZ,EZ) ## to be double checked !!
    FE = doubleCheckFaceBoundaries(FE,Z,FZ,EZ)
    
    # remove 0 indices from FE relation
    FE = [[f  if f!=0 else 1 for f in face] for face in FE]
    EF_angle = faceSlopeOrdering(model,FE,Z)
    
    V,CV,FV,EV,CF,CE,COE = facesFromComponents((Z,FZ,EZ),FE,EF_angle)
    return V,CV,FV,EV,CF,CE,COE,FE

