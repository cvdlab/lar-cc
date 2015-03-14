""" Module for Boolean computations between geometric objects """
from pyplasm import *
""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')
from inters import *
DEBUG = True

""" Coding utilities """
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
    transformMat = mat( newFacet[1:d] + [covector[1:]] ).T.I
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
        print "k,edges =",k,edges
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
        print "\ne,edgeFaces,EW[e] =", e,edgeFaces,EW[e]
        for f in edgeFaces:
            print "f,FW[f],TW[f] =", f,FW[f],TW[f]
            for t in TW[f]:
                if set(EW[e]).intersection(t)==set(EW[e]):
                    ET[e] += [t]
                    print "ET[e] =",ET[e]
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
from bool1 import invertRelation

def faceSlopeOrdering(model):
    V,FV,EV = model

    print "\n V,FV,EV =",V,FV,EV

    FE = crossRelation(FV,EV)

    print "\n FE =",FE

    EF,EF_angle = invertRelation(FE),[]

    print "\n EF,EF_angle =",EF,EF_angle

    for e,ef in enumerate(EF):

        print "\n e,ef =",e,ef

        v1,v2 = EV[e]

        print "\n v1,v2 =",v1,v2

        E = array(UNITVECT(VECTDIFF([V[v2],V[v1]])))

        print "\n E =",E

        refFace = EF[e][0]

        print "\n refFace =",refFace

        ef_angle = []

        print "\n ef_angle =",ef_angle

        for face in ef:

            print "\n face =",face

            fverts = [ VECTDIFF([ V[k],V[v1] ]) for k in FV[face] if k!=v1 ]

            print "\n fverts =",fverts

            fvects = array([ E.dot( vf ) for vf in fverts ])

            print "\n fvects =",fvects

            index_min = fvects.argmin()

            print "\n index_min =",index_min

            W = fverts[index_min]

            print "\n W =",W

            if not verySmall(E.dot(W)): W = UNITVECT(orthoProject(list(E))(W))

            print "\n W =",W

            if (face == refFace): refVect = array(UNITVECT(W))

            print "\n refVect =",refVect

            angle = math.atan2(VECTNORM(VECTPROD([refVect,W])),refVect.dot(W))

            print "\n angle =",angle

            ef_angle += [180*angle/PI]

            print "\n ef_angle =",ef_angle

        pairs = sorted(zip(ef_angle,ef))

        print "\n pairs =",pairs

        EF_angle += [[pair[1] for pair in pairs]]

        print "\n EF_angle =",EF_angle

    return EF_angle

