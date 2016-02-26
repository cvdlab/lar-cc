""" Module for Boolean computations between geometric objects """
from larlib import *
import inters,triangulation
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
import larcc

def MKTRIANGLES(model,color=False):
    V,FV,EV = model
    lenV = len(V)
    VV = AA(LIST)(range(len(V)))
    FE = larcc.crossRelation(FV,EV,VV)
    if len(V[0]) == 2: V=[v+[0] for v in V]
    triangleSets = boundaryTriangulation(V,FV,EV,FE)
    if color:
        colors = [CYAN,MAGENTA,WHITE,RED,YELLOW,GREEN,GRAY,ORANGE,BLACK,BLUE,PURPLE,BROWN]
        return [ COLOR(colors[k%12])(STRUCT([MKPOL([verts,[[3,2,1]],None]) 
            for verts in triangledFace])) for k,triangledFace in enumerate(triangleSets) ]
    else:
        return [ STRUCT([MKPOL([verts,[[3,2,1]],None]) for verts in triangledFace])
                for triangledFace in triangleSets ]

def MKSOLID(*model): 
    V,FV,EV = model
    VV = AA(LIST)(range(len(V)))
    FE = larcc.crossRelation(FV,EV,VV)
    pivot = V[0]
    #VF = invertRelation(FV) 
    #faces = [face for face in FV if face not in VF[0]]
    faces = [face for face in FV]
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
        below,above = inters.splitOnThreshold(boxes,bucket,1)
        below1,above1 = inters.splitOnThreshold(boxes,above,2)
        below2,above2 = inters.splitOnThreshold(boxes,below,2) 
               
        below11,above11 = inters.splitOnThreshold(boxes,above1,3)
        below21,above21 = inters.splitOnThreshold(boxes,below1,3)        
        below12,above12 = inters.splitOnThreshold(boxes,above2,3)
        below22,above22 = inters.splitOnThreshold(boxes,below2,3)  
              
        inters.splitting(above1,below11,above11, finalBuckets,splittingStack)
        inters.splitting(below1,below21,above21, finalBuckets,splittingStack)
        inters.splitting(above2,below12,above12, finalBuckets,splittingStack)
        inters.splitting(below2,below22,above22, finalBuckets,splittingStack)
        
        finalBuckets = list(set(AA(tuple)(finalBuckets)))
    parts = inters.geomPartitionate(boxes,finalBuckets)
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
import larcc
def submodel(V,FV,EV):
    lenV = len(V)
    VV = AA(LIST)(range(lenV))
    FE = larcc.crossRelation0(lenV,FV,EV)
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
def veryClose(p,q):
    out = False
    if VECTNORM(VECTDIFF([p,q])) <= 0.002:
        out = True
    return out

def removeExternals(M,V,EV,fe, z,fz,ez):
    (Z,FZ,EZ) = copy((z,fz,ez))
    inputFace = Struct([(V,[EV[e] for e in fe])])
    verts,ev = larApply(M)(struct2lar(inputFace))
    pol = [eval(vcode(vert[:-1])) for vert in verts],ev
    out = []
    classify = triangulation.pointInPolygonClassification(pol)
    for k,point in enumerate(Z):
        if classify(point)=="p_out":  out += [k]

    # verify all v in out w.r.t. pol[0]
    trueOut = []
    for v in out: 
        onBoundary = False
        for p in pol[0]:
            if veryClose(Z[v],p):
                print "v,p,'close'", Z[v],p
                onBoundary = True
                Z[v] = p
        if not onBoundary: trueOut += [v]
    
    FW = [f for f in FZ if not any([v in trueOut for v in f])]  # trueOut
    EW = [e for e in EZ if not any([v in trueOut for v in e])]  # trueOut
    return Z,FW,EW

""" Space partitioning via submanifold mapping """
import larcc
from larcc import *

def spacePartition(V,FV,EV, parts):
    VV = AA(LIST)(range(len(V)))
    FE = larcc.crossRelation(FV,EV,VV)
    submodel0 = submodel(V,FV,EV)
    out = []
    
    """ input: face index f; candidate incident faces F[f]; """
    for f,F in enumerate(parts):
        """ Selection of the LAR submodel S(f) := (V,FV,EV)(f) restricted to [f]+F[f] """    
        fF,fE = submodel0(f,F)
        subModel = Struct([(V,[FV[g] for g in fF],[EV[h] for h in fE])])
        sV,sFV,sEV = struct2lar(subModel)
        #VIEW(STRUCT(MKPOLS((sV,sFV+sEV)) + [COLOR(RED)(STRUCT(MKPOLS((V,[FV[f]]))))]))
        
        """ Computation of submanifold map M moving f to z=0 """
        pivotFace = [V[v] for v in FV[f]]
        M = submanifoldMapping(pivotFace)  # submanifold transformation
        
        """ Transformation of S(f) by M, giving S = (sW,sEW) := M(S(f)) """
        sW,sFW,sEW = larApply(M)((sV,sFV,sEV))
        
        red = COLOR(RED)(STRUCT(MKPOLS(larApply(M)((V,[FV[f]])))))
        #VIEW(STRUCT(MKPOLS((sW,sFW+sEW)) + [red]))
        
        """ filtering of EW edges traversing z=0, giving EZ edges and incident faces FZEZ """
        sFE = larcc.crossRelation0(len(V),sFW,sEW)    
        edges = list(set([ e for k,face in enumerate(sFW)  for e in sFE[k] 
                    if meetZero(sW, sEW[e]) ]))
        edgesPerFace = [ [e for e in sFE[k] if meetZero(sW, sEW[e])] 
                    for k,face in enumerate(sFW) ]
        edges = list(set(CAT(edgesPerFace)))
        WW = AA(LIST)(range(len(sW)))
        wires = [sEW[e] for e in edges]
        #VIEW(STRUCT(MKPOLS((sW,wires)) + [red]))
        
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
        lines = [[line[0],line[-1]]  if len(line)>2 else line for line in lines]
        #VIEW(STRUCT(AA(POLYLINE)(lines) + [red]))
                
        """ Construction of the planar set FX,EX of faces and lines """
        u,fu,eu,polygons = inters.larFromLines(lines)
        #VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((u,fu+eu))))
        
        """ Remove external vertices """
        w,fw,ew = removeExternals(M,V,EV,FE[f], u,fu,eu)
        #VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((w,fw+ew))))
        struct = Struct([(w,fw,ew)])
        z,fz,ez = struct2lar(struct)
        
        #VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((z,fz+ez))))
        
        w,fw,ew = larApply(M.I)(([v+[0.0] for v in z],fz,ez))
        out += [Struct([(w,fw,ew)])]
    return struct2lar(Struct(out))

""" 3D boundary triangulation of the space partition """

def orientTriangle(pointTriple):
    v1 = array(pointTriple[1])-pointTriple[0]
    v2 = array(pointTriple[2])-pointTriple[0]
    if cross(v1,v2)[2] < 0: return REVERSE(pointTriple)
    else: return pointTriple
    
from copy import copy
from scipy import spatial

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
              
        # Construction of CDT (Constrained Delaunay Triangulation) for LAR face
        model = (verts2D,[range(len(verts2D))],EW)
        struct = Struct([model])
        U,EU = triangulation.structBoundaryModel(struct)
        W,FW,EW,polygons = inters.larFromLines([[U[u],U[v]] for u,v in EU])

        def filtering(polygons):
            if len(polygons)>1 and max(AA(len)(polygons))==1:
                poly = sorted(CAT(polygons),key=len,reverse=True)
                poly = [poly[0] for comp in poly[1:] if set(comp).issubset(set(poly[0]))]
                return [poly]
            else: return polygons

        thePolygons = filtering(polygons)
        setsOfTriangles = triangulation.polygons2TriangleSet(W,thePolygons)
        trias = [[p+[1],q+[1],r+[1]] for p,q,r in setsOfTriangles[0]]
        
        inverseMap = transform.I
        trias = AA(mapVerts(inverseMap))(trias)
        triangledFace += [[v[:-1] for v in triangle] for triangle in trias]
        triangleSet += [triangledFace]
    return triangleSet

def triangleIndices(triangleSet,W):
    tree = spatial.cKDTree(W)
    TV,FT,t = [],[],-1
    for face in triangleSet:
        ft = []
        for triangle in face:
            vertices = tree.query(triangle,1)[1].tolist()
            t += 1
            TV += [vertices]
            ft += [t]
        FT += [ft]
    VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,TV))))
    return TV,FT


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

def faceSlopeOrdering(model,FE):
    V,FV,EV = model
    VIEW(COLOR(YELLOW)(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,EV)))))
    triangleSet = boundaryTriangulation(V,FV,EV,FE) # corrected with non-contractible faces
    VIEW(EXPLODE(1.2,1.2,1.2)(AA(JOIN)( AA(POLYLINE)(CAT(triangleSet)) )))
    VIEW(EXPLODE(1.2,1.2,1.2)( AA(POLYLINE)(AA(lambda tri: tri+[tri[0]])(CAT(triangleSet))) ))
    TV,FT = triangleIndices(triangleSet,V) 
    VV = AA(LIST)(range(len(V)))
    TE = crossRelation(TV,EV,VV)
    ET,ET_angle = invertRelation(TE),[]
    #import pdb; pdb.set_trace()
    for e,et in enumerate(ET):
        v1,v2 = EV[e]
        v1v2 = set([v1,v2])
        et_angle = []
        t0 = et[0]
        tverts = [v1,v2] + list(set(TV[t0]).difference(v1v2))
        e3 = UNITVECT(VECTDIFF([ V[tverts[1]], V[tverts[0]] ]))
        e1 = UNITVECT(VECTDIFF([ V[tverts[2]], V[tverts[0]] ]))
        e2 = cross(array(e1),e3).tolist()
        basis = mat([e1,e2,e3]).T
        transform = basis.I
        normals = []
        Tvs = []
        for triangle in et:
            verts = TV[triangle]
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
    EF_angle = ET_to_EF_incidence(TV,FV,FT, ET_angle)
    return EF_angle, ET,TV,FT

""" Edge-triangles to Edge-faces incidence """
def ET_to_EF_incidence(TW,FW,FT, ET_angle):
    tableFT = FT
    tableTF = invertRelation(tableFT)
    EF_angle = [[tableTF[t][0] for t in triangles] for triangles in ET_angle]
    #assert( len(EF_angle) == 2*len(FW) )
    return EF_angle




""" Cells from $(d-1)$-dimensional LAR model """
def facesFromComponents(model,FE,EF_angle):

    accumulated = []
    def viewStep (CF,CV,CE,COE,accumulated):
        VV = AA(LIST)(range(len(V)))
        edges = list(set(CE[-1]).difference(accumulated))
        accumulated = CE[-1]
        submodel = STRUCT(MKPOLS((V,[EV[k] for k in edges])))
        VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,1))

    print "\nECCOMI\n"
    # initialization
    V,FV,EV = model
    visitedCell = [[ None, None ] for k in range(len(FV)) ]
    face = 0
    boundaryLoop,_ = boundaryCycles(FE[face],EV)
    boundaryLoop = boundaryLoop[0]
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
    viewStep (CF,CV,CE,COE,accumulated)
    
    # main loop
    #import pdb; pdb.set_trace()
    while True:
        face, edge = startCell(visitedCell,FE,EV)
        print "face, edge =",face, edge
        if face == -1: break
        boundaryLoop,_ = boundaryCycles(FE[face],EV)
        boundaryLoop = boundaryLoop[0]
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
        viewStep (CF,CV,CE,COE,accumulated)
    return V,CV,FV,EV,CF,CE,COE

""" Cycles orientation """
def cyclesOrient(pcycles,fcycle,EV):
    if set(AA(ABS)(pcycles)).difference(fcycle)==set(): return []
    ofcycle,_ = boundaryCycles(fcycle,EV) # oriented 
    ofcycle = ofcycle[0] # oriented 
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
    faceLoop,_ = boundaryCycles(FE[face],EV)
    faceLoop = faceLoop[0]
    commonEdges = set(faceLoop).intersection(boundaryLoop)
    if commonEdges == set() or commonEdges == {0}: 
        faceLoop = reverseOrientation(faceLoop)
        commonEdges = set(faceLoop).intersection(boundaryLoop)
    theEdge = list(commonEdges)[0]
    #if theEdge==0: theEdge = list(commonEdges)[1]
    return -theEdge,face

""" Extend LAR edges with a given (2D) triangulation """
from collections import OrderedDict

def extendEV(EV,ET,TV):
    EVdict = OrderedDict([(edge,k) for k,edge in enumerate(EV)])
    n = len(EV)-1
    for e,(u,w) in enumerate(EV):
        for t in ET[e]:
            v1,v2,v3 = TV[t]
            v = list({v1,v2,v3}.difference([u,w]))[0]
            if u<v: newEdge = (u,v)
            
            else: newEdge = (v,u)
            if not newEdge in EVdict: 
                n += 1
                EVdict[newEdge]=n
                
            if w<v: newEdge = (w,v)
            else: newEdge = (v,w)
            if not newEdge in EVdict: 
                n += 1
                EVdict[newEdge]=n
    return EVdict.keys()

""" Signed boundary operator for a general LAR 2-complex """
def larSignedBoundary(larModel,triaModel,FT):
    inputOp = signedSimplicialBoundary(*triaModel)
    outputOp = boundary(*larModel)
    (n,m),(p,q) = inputOp.shape, outputOp.shape
    
    for h in range(p):   # for each LAR face
        for k,triangles in enumerate(FT):
            val = [inputOp[h,t] for t in triangles if inputOp[h,t]!=0.0]
            if val!=[]: outputOp[h,k] = val[0]
    return outputOp

""" Get single solid cell """
def getSolidCell(FE,face,visitedCell,boundaryLoop,EV,EF_angle,V,FV):

    def orientFace(face,boundaryLoop): 
        for e in boundaryLoop:
            if ABS(e) in FE[face]: return -e

    coe = [orientFace(face,boundaryLoop)]
    cf = [face] 
    #import pdb; pdb.set_trace()
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
        VIEW(STRUCT( MKPOLS((V,[EV[h] for f in cf for h in FE[f]])) )) #add EV!
    return cf,coe

""" Main procedure of arrangement partitioning """
import inters

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
    Z,FZ,EZ = inters.larSimplify((Z,FZ,EZ),radius=0.0001)
    model = Z,FZ,EZ

    ZZ = AA(LIST)(range(len(Z)))
    submodel = STRUCT(MKPOLS((Z,EZ)))
    VIEW(larModelNumbering(1,1,1)(Z,[ZZ,EZ,FZ],submodel,0.6)) 

    ZZ = AA(LIST)(range(len(Z)))
    FE = crossRelation(FZ,EZ,ZZ) ## to be double checked !!
    EF_angle, ET,TV,FT = faceSlopeOrdering(model,FE)
    
    V,CV,FV,EV,CF,CE,COE = cellsFromComponents((Z,FZ,EZ),FE,EF_angle, ET,TV,FT)
    V,CV,FV,EV,CF,CE,COE = facesFromComponents((Z,FZ,EZ),FE,EF_angle)
    return V,CV,FV,EV,CF,CE,COE,FE

