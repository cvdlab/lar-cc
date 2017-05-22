""" Module for Boolean computations between geometric objects """
from larlib import *
import inters,triangulation
from copy import copy
DEBUG = False

""" Coding utilities """
global count
""" Generation of a random 3D point """
def rpoint3d():
    return eval( vcode(4)([ random.random(), random.random(), random.random() ]) )

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
    return tuple(eval(vcode(4)(v1))), tuple(eval(vcode(4)(v2))), tuple(eval(vcode(4)(v3)))

""" Containment boxes """
def containmentBoxes(randomPointArray,qualifier=0):
    if len(randomPointArray[0])==2:
        boxes = [eval(vcode(4)([min(x1,x2), min(y1,y2), min(z1,z2), 
                             max(x1,x2), max(y1,y2), max(z1,z2)]))+[qualifier]
                for ((x1,y1,z1),(x2,y2,z2)) in randomPointArray]
    elif len(randomPointArray[0])==3:
        boxes = [eval(vcode(4)([min(x1,x2,x3), min(y1,y2,y3), min(z1,z2,z3), 
                             max(x1,x2,x3), max(y1,y2,y3), max(z1,z2,z3)]))+[qualifier]
                for ((x1,y1,z1),(x2,y2,z2),(x3,y3,z3)) in randomPointArray]
    elif len(randomPointArray[0])==4:
        boxes = [eval(vcode(4)([min(x1,x2,x3,x4), min(y1,y2,y3,y4), min(z1,z2,z3,z4), 
                             max(x1,x2,x3,x4), max(y1,y2,y3,y4), max(z1,z2,z3,z4)]))+[qualifier]
                for ((x1,y1,z1),(x2,y2,z2),(x3,y3,z3),(x4,y4,z4)) in randomPointArray]
    return boxes

def containmentBoxes(randomPointArray,qualifier=0):
    def minmax(pointArray):
        coords = TRANS(pointArray)
        return AA(min)(coords) + AA(max)(coords)
    boxes = [eval(vcode(4)( minmax(pointArray) ))+[qualifier] for pointArray in randomPointArray]
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
        boxes += [eval(vcode(4)([min(x1,x2),min(y1,y2),min(z1,z2),max(x1,x2),max(y1,y2),max(z1,z2)]))+[(qualifier,k)]]
    return boxes

""" Generation of a list of HPCs from a LAR model with non-convex faces """
import boundary

def MKTRIANGLES(model,color=False):
    V,FV,EV = model
    lenV = len(V)
    VV = AA(LIST)(range(len(V)))
    
    boundaryOperator = boundary.larSignedBoundary2(V,FV,EV)
    FEbasis = boundary.signedBasis(boundaryOperator)
    FE,signs = TRANS(FEbasis)
    
    if len(V[0]) == 2: V=[v+[0] for v in V]
    triangleSets = boundaryTriangulation(V,FV,EV,FE)
    triangleSet = [[(p1,p2,p3) if sign==1 else (p2,p1,p3) for p1,p2,p3 in triangleSet] 
                    for sign,triangleSet in zip(signs,triangleSets)]
    if color:
        colors = [CYAN,MAGENTA,WHITE,RED,YELLOW,GREEN,GRAY,ORANGE,BLACK,BLUE,PURPLE,BROWN]
        return [ COLOR(colors[k%12])(STRUCT([MKPOL([verts,[[1,2,3]],None]) 
            for verts in triangledFace])) for k,triangledFace in enumerate(triangleSets) ]
    else:
        return [ STRUCT([MKPOL([verts,[[1,2,3]],None]) for verts in triangledFace])
                for triangledFace in triangleSets ]
"""
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
"""

""" Utility to transform a dictionary to a function on the keys """
def dict2fun(d):
    def dict2fun0(k): return d[k]
    return dict2fun0 

""" Coherent orientation of boundary 2-faces """
def boundaryOrientation(V,EV,FEbasis,signedBoundary,triangleSets):            
    vdict = OrderedDict([ (vcode(3)(v),k) for k,v in enumerate(V)  ])
    edict = OrderedDict([ (edge,k) for k,edge in enumerate(EV) ])
    triaVerts = [[AA(dict2fun(vdict))(AA(vcode(3))(t)) for t in triangleSet] 
                 for triangleSet in triangleSets]
    triaEdges = []
    for triangles in triaVerts:
        tria2edgs = []
        for v0,v1,v2 in triangles:
            tris = []
            if (v0,v1) in edict: tris += [edict[(v0,v1)]]
            if (v1,v2) in edict: tris += [edict[(v1,v2)]]
            if (v2,v0) in edict: tris += [edict[(v2,v0)]]
            tria2edgs += [tris]
        triaEdges += [CAT(tria2edgs)[0] if CAT(tria2edgs)[0]!=0 else CAT(tria2edgs)[1]]
    
    orientations = zip(triaEdges,[[sign,FEbasis[f]] for f,sign in signedBoundary])
    
    tests = [set([triaEdges]).intersection((sign*array(signs)*array(edgechain)).tolist()) 
        for triaEdges,(sign,(edgechain,signs)) in orientations]
        
    theSigns = [1 if val!= set([]) else -1 for val in tests]
    return triangleSets,theSigns

""" Generation of LAR B-rep from a LAR model with non-convex faces """

def BREP (model,color=False):
    # intrinsic orientation of input 2-faces
    V,FV,EV = model
    VV = AA(LIST)(range(len(V)))
    boundaryOperator = boundary.larSignedBoundary2(V,FV,EV)
    FEbasis = boundary.signedBasis(boundaryOperator)
    FE,facesigns = TRANS(FEbasis)
    
    # computation of boundary 2-faces
    csrBoundary3,CF,_ = boundary.larSignedBoundary3((V,FV,EV))
    
    signedBoundary = csrBoundary3 * (len(CF)*[[1]])  ## TODO !!
    
    cells,_ = TRANS(signedBoundary)
    fv = [FV[f] for f in cells]
    ev = [EV[e] for e in set(CAT([[e for e in FE[f]] for f in cells])) ]
    fe = larcc.crossRelation(fv,ev,VV)
    triangleSets = boundaryTriangulation(V,fv,ev,fe) 
    
    # computation of coherent orientation of boundary 2-faces
    triangleSets,theSigns = boundaryOrientation(V,EV,FEbasis,signedBoundary,triangleSets)

    # visualization of boundary
    if color:
        colors = [CYAN,MAGENTA,WHITE,RED,YELLOW,GREEN,GRAY,ORANGE,BLACK,BLUE,PURPLE,BROWN]
        return [ COLOR(colors[k%12])(
                 STRUCT([MKPOL([verts,[[3,2,1]],None]) for verts in triangledFace])
             if sign==-1 else 
             STRUCT([MKPOL([verts,[[1,2,3]],None])  for verts in triangledFace])
                 ) 
        for k,(sign,triangledFace) in enumerate(zip(theSigns,triangleSets)) ]
    else:
        return [ STRUCT([MKPOL([verts,[[3,2,1]],None])  for verts in triangledFace])
             if sign==-1 else 
             STRUCT([MKPOL([verts,[[1,2,3]],None])  for verts in triangledFace])
        for sign,triangledFace in zip(theSigns,triangleSets) ]
""" Generation of LAR B-rep from a LAR model with non-convex faces """
if __name__=="__main__":

    V,[VV,EV,FV,CV] = larCuboids([2,2,2],True)
    cubeGrid = Struct([(V,FV,EV)],"cubeGrid")
    cubeGrids = Struct(2*[cubeGrid,t(.5,.5,.5),r(0,0,PI/6)])

    V,FV,EV = struct2Marshal(cubeGrids)
    VIEW(EXPLODE(1.2,1.2,1.2)(BREP((V,FV,EV),color=True) ))
    VIEW(EXPLODE(1.2,1.2,1.2)(BREP((V,FV,EV)) ))
    VIEW(STRUCT(BREP((V,FV,EV)) ))



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
    #print "\nfacet =",facet
    covector = COVECTOR(facet)
    #print "\ncovector =",covector
    translVector = facet[0]
    #print "translVector =",translVector
    # translation 
    newFacet = [ VECTDIFF([v,translVector]) for v in facet ]
    #print "newFacet =",newFacet
    # linear transformation: boundaryFacet -> standard (d-1)-simplex
    if isclose(0.,covector[1]) and isclose(0.,covector[2]): ## x and y components ! (hpc format)
        m = mat(np.eye(3))
    else:
        d = len(facet[0])
        #print "d =",d
        m = mat( newFacet[1:d] + [covector[1:]] )
        #print "m =",m
        if abs(det(m))<0.0001:
            for k in range(len(facet)-2):
                m = mat( newFacet[1+k+1:d+k+1] + [covector[1:]] )
                #print "\nm =",m
                if abs(det(m))>0.0001: break
    transformMat = m.T.I
    #print "transformMat =",transformMat
    # transformation in the subspace x_d = 0
    out = (transformMat * (mat(newFacet).T)).T.tolist()
    #print "out =",out
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

def submanifoldMapping(pivotFace):
    tx,ty,tz = pivotFace[0]
    transl = mat([[1,0,0,-tx],[0,1,0,-ty],[0,0,1,-tz],[0,0,0,1]])
    facet = [ VECTDIFF([v,pivotFace[0]]) for v in pivotFace ]
    normal = UNITVECT(COVECTOR(facet)[1:])
    a = normal
    b = [0,0,1]
    axis = UNITVECT(VECTPROD([a,b]))
    angle = math.atan2(VECTNORM(cross(a,b)), dot(a,b))    
    
    # general 3D rotation (Rodrigues' rotation formula)    
    m = scipy.identity(4)
    cos = COS(angle); sin = SIN(angle)
    I = scipy.identity(3) ; u = axis
    Ux = scipy.array([
        [0,        -u[2],      u[1]],
        [u[2],        0,     -u[0]],
        [-u[1],     u[0],         0]])
    UU = scipy.array([
        [u[0]*u[0],    u[0]*u[1],    u[0]*u[2]],
        [u[1]*u[0],    u[1]*u[1],    u[1]*u[2]],
        [u[2]*u[0],    u[2]*u[1],    u[2]*u[2]]])
    m[:3,:3] = cos*I + sin*Ux + (1.0-cos)*UU
    
    mapping = mat(m)
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
import boundary,inters

def veryClose(edge,p):
    ((x1,y1),(x2,y2)),(x0,y0) = edge,p
    distance = ABS((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1) / VECTNORM(VECTDIFF([(x1,y1),(x2,y2)]))
    if distance <= 0.00001: return True
    return False
    
def removeExternals(M,V,EV,fe,fv, z,fz,ez):
    w,fw,ew = struct2lar(Struct([larApply(M)((V,[fv],[EV[e] for e in fe]))])) # part mapped to 2D
    newEdges = boundary.larOffset2D(([v[:-1] for v in w],fw,ew),offset=0.0001)

    #w,fw,ew,_ = inters.larFromLines(newEdges)
    w,ew = bruteForceIntersect(newEdges)
    w,polygons,ew = triangulation.larPair2Triple((w,ew))
    fw = AA(list)(AA(set)(AA(CAT)(polygons)))

    pol = w,ew
    out = []
    classify = triangulation.pointInPolygonClassification(pol)
    for k,point in enumerate(z):
        if classify(point)=="p_out":  out += [k]

    # verify all v in out w.r.t. pol[0]
    trueOut = []
    w = pol[0]
    for k in out: 
        p = z[k]
        onBoundary = False
        for u,v in pol[1]:
            if veryClose((w[u],w[v]),p):
                onBoundary = True
                z[k] = p
        if not onBoundary: trueOut += [k]
    
    fw = [f for f in fz if not any([v in trueOut for v in f])]  # trueOut
    ew = [e for e in ez if not any([v in trueOut for v in e])]  # trueOut
    return z,fw,ew
        
#VIEW(STRUCT(AA(MK)([z[k] for k in set(range(len(z))).difference(out)])))
#VIEW(STRUCT(AA(MK)(z)))
#VIEW(EXPLODE(1.2,1.2,1)(AA(COLOR(CYAN))(MKPOLS((z,ez)))+AA(COLOR(YELLOW))(MKPOLS((w,ew)))))

""" Compute face intersections with z=0 """
def computeCrossingLines(edges,sW,sFW,sEW,sFE):
   crossEdges = []
   def isClose(a,b): return abs(a-b)<10**-5
   for e in edges:
      z1,z2 = sW[sEW[e][0]][2], sW[sEW[e][1]][2]
      cross = False
      if isClose(z1,z2):   cross=False
      elif z1*z2 < 0:      cross=True
      elif isClose(0.,z1*z2): cross=True
      else: cross = False
      if cross: crossEdges += [e]   
   if crossEdges != []:
      #VIEW(STRUCT(MKPOLS((sW,[sEW[e] for e in crossEdges]))))
      sEF = invertRelation(sFE)
      crossFaces = list(set(CAT([sEF[e] for e in crossEdges])))
      #VIEW(STRUCT(MKPOLS((sW,[sFW[f] for f in crossFaces]+[sEW[e] for e in crossEdges]))))
      #VIEW(STRUCT(MKPOLS((sW,[sEW[e] for e in crossEdges]))))
      edgeCrossSets = [list(set(crossEdges).intersection(sFE[f])) for f in crossFaces]    
      
      def points2lines(pointSet):
         #  preconditions:
         #  1. len(pointSet) > 2;
         #  2. 2D points in pointSet are aligned.
         pointSet = list(set(AA(tuple)(AA(eval)(AA(vcode(6))(pointSet)))))
         # postcondition: pointSet has even length
         if len(pointSet) == 2: return [pointSet]
         lines = []
         p1,p2 = pointSet[:2]
         # TODO: multiple lines
         return lines

      pointSets2d = []
      for Set in edgeCrossSets:
         pts2d = []
         for edge in Set:
            v1,v2 = [sW[v] for v in sEW[edge]]
            [x1,y1,z1], [x2,y2,z2] = v1, v2
            u = z1/(z1-z2)
            p2d = x1 + u*(x2-x1), y1 + u*(y2-y1)
            pts2d += [p2d]
         pointSets2d += [pts2d]

      crossingLines = []
      for pointSet in pointSets2d:
         if len(pointSet)==2:
            crossingLines += [pointSet]
         else:
            crossingLines += points2lines(pointSet)

      return crossingLines
   else: return []

""" Brute-force intersection of 2D lines """
def bruteForceIntersect(lines):
   n = len(lines)
   #transform data
   lines = [[eval(vcode(4)(p)) for p in line] for line in lines]
   #end transform
   verts = list(set([tuple(v) for line in lines for v in line]))
   vertdict = OrderedDict([(key,k) for k,key in enumerate(verts)])
   EV = [[vertdict[tuple(p)] for p in line] for line in lines]
   pairs = [(h,k) for h in range(n) for k in range(h+1,n)]
   linepairs = [[lines[h],lines[k]] for h,k in pairs]
   # prepare data for line pairs
   linedata = [[ax,ay,bx,by,cx,cy,dx,dy] 
      for [[(ax,ay),(bx,by)],[(cx,cy),(dx,dy)]] in linepairs]
   # assemble intersection determinants
   determinants = [ det(mat([[ax-bx,dx-cx], [ay-by,dy-cy]])) 
      for [ax,ay,bx,by,cx,cy,dx,dy] in linedata]
   # parameter pairs by Cramer's rule (for oriented edges of f face)
   alpha = [det(mat([[dx-bx,dx-cx],[dy-by,dy-cy]]))/D  if abs(D)>.00001 else 0 
      for D,(ax,ay,bx,by,cx,cy,dx,dy) in zip(determinants,linedata)]
   beta = [det(mat([[ax-bx,dx-bx],[ay-by,dy-by]]))/D  if abs(D)>.00001 else 0 
      for D,(ax,ay,bx,by,cx,cy,dx,dy) in zip(determinants,linedata)]
   # intersection points
   vdata,edata = defaultdict(list),defaultdict(list)
   newverts = [ tuple(AA(COMP([tuple,eval,vcode(4)]))([ 
      (a*mat(p1)+(1-a)*mat(p2)).tolist()[0], 
      [a,b,h,k] ]))
      for (a,b,(h,k)),[[p1,p2],[q1,q2]] in zip(zip(alpha,beta,pairs),linepairs) 
      if 0<=a<=1 and 0<=b<=1 ]
   for vert,datum in newverts:
      vdata[vert] += [datum]
   for k,(key,datum) in enumerate(vdata.items()):
      for a,b,edge1,edge2 in datum:
         edata[int(edge1)] += [a]
         edata[int(edge2)] += [b]
   edgeParameters = [sorted(set(edge)) for k,edge in edata.items()]
   points = [[(a*mat(verts[p])+(1-a)*mat(verts[q])).tolist()[0] 
         for a in params] for params,(p,q) in zip(edgeParameters,EV)]
   m = len(vertdict)
   for point in CAT(points):
      vertex = tuple(eval(vcode(4)(point)))
      if not vertex in vertdict:
         vertdict[vertex] = m
         m += 1
   edgeVerts = [[vertdict[tuple(eval(vcode(4)(point)))] for point in edge] 
      for edge in points]
   V = AA(list)(vertdict.keys())
   edges = [[[v,part[k+1]] for k,v in enumerate(part[:-1])] for part in edgeVerts]
   EV = sorted(set(AA(tuple)(AA(sorted)(CAT(edges)))))
   return V,EV

""" Space partitioning via submanifold mapping """
import larcc,inters,triangulation
from larcc import *

def lineExtend(epsilon):
   def lineExtend1(line):
      ((x1,y1),(x2,y2)) = line
      def x(a): return x1 + a*(x2-x1)
      def y(a): return y1 + a*(y2-y1)
      return ((x(-epsilon),y(-epsilon)),(x(1+epsilon),y(1+epsilon)))
   return lineExtend1

def spacePartition(V,FV,EV, parts):
    VV = AA(LIST)(range(len(V)))
    FE = larcc.crossRelation(FV,EV,VV)
    submodel0 = submodel(V,FV,EV)
    out = []
    def isClose(a,b): return abs(a-b)<10**-5
    
    """ input: face index f; candidate incident faces F[f]; """
    for f,F in enumerate(parts):
        print "\n\nf,F =",f,F
        """ Selection of the LAR submodel S(f) := (V,FV,EV)(f) restricted to [f]+F[f] """    
        fF,fE = submodel0(f,F)
        sbModel = Struct([(V,[FV[g] for g in fF],[EV[h] for h in fE])])
        sV,sFV,sEV = struct2lar(sbModel)
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
        alledges = list(set([ e for k,face in enumerate(sFW)  for e in sFE[k] 
                    if meetZero(sW, sEW[e]) ]))
        edgesPerFace = [ [e for e in sFE[k] if meetZero(sW, sEW[e])] 
                    for k,face in enumerate(sFW) ]
        edges = list(set(CAT(edgesPerFace)))
        
        wires = [sEW[e] for e in edges]
        #VIEW(STRUCT(MKPOLS((sW,wires)) + [red]))
        
        """ for each face in FZEZ, computation of the aligned set of points p(z=0) """
        
        edges2D = [[w1,w2] for w1,w2 in wires if (isClose(0,sW[w1][2]) and isClose(0,sW[w2][2]))]
        lines2d = [[sW[u][:-1],sW[v][:-1]] for u,v in edges2D] + computeCrossingLines(
                    edges,sW,sFW,sEW,sFE)
                    
        #lines2D = AA(lineExtend(0.0001))(lines2d)
        #lines2D = AA(AA(COMP([eval,vcode(0.0000000000001)])))(lines2d)          

        #VIEW(STRUCT(AA(POLYLINE)(lines2D) ))#+ [red]))
        
        u,fu,eu,_ = inters.larFromLines(lines2d)
        u,polygons,eu = triangulation.larPair2Triple((u,eu))
        fu = AA(list)(AA(set)(AA(CAT)(polygons)))
        
        z,fz,ez = u,fu,eu
        
        #Remove external vertices 
        #z,fz,ez = removeExternals(M,V,EV,FE[f],FV[f], u,fu,eu)  # BUG !!!!  <<<<<<<<<<
        #VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((z,fz+ez))))
        
        
        w,fw,ew = larApply(M.I)(([v+[0.0] for v in z],fz,ez))
        #VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((w,fw+ew))))
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
        if setsOfTriangles != []:
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
    #VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,TV))))
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
    #VIEW(COLOR(YELLOW)(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,EV)))))
    triangleSet = boundaryTriangulation(V,FV,EV,FE) # corrected with non-contractible faces
    #VIEW(EXPLODE(1.2,1.2,1.2)(AA(JOIN)( AA(POLYLINE)(CAT(triangleSet)) )))
    #VIEW(EXPLODE(1.2,1.2,1.2)( AA(POLYLINE)(AA(lambda tri: tri+[tri[0]])(CAT(triangleSet))) ))
    TV,FT = triangleIndices(triangleSet,V) 
    VV = AA(LIST)(range(len(V)))
    TE = crossRelation(TV,EV,VV)
    ET,ET_angle = invertRelation(TE),[]
    #import pdb; pdb.set_trace()
    for e,et in enumerate(ET):
        if et != []:
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
    import triangulation

    accumulated = []
    def viewStep (CF,CV,CE,COE,accumulated):
        VV = AA(LIST)(range(len(V)))
        edges = list(set(CE[-1]).difference(accumulated))
        accumulated = CE[-1]
        submodel = STRUCT(MKPOLS((V,[EV[k] for k in edges])))
        #VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,1))

    # initialization
    V,FV,EV = model
    visitedCell = [[ None, None ] for k in range(len(FV)) ]
    face = 0
    boundaryLoop,_ = triangulation.boundaryCycles(FE[face],EV)
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
        if face == -1: break
        boundaryLoop,_ = triangulation.boundaryCycles(FE[face],EV)
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
import triangulation
def cyclesOrient(pcycles,fcycle,EV):
    if set(AA(ABS)(pcycles)).difference(fcycle)==set(): return []
    ofcycle,_ = triangulation.boundaryCycles(fcycle,EV) # oriented 
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
import triangulation

def reverseOrientation(chain):
    return REVERSE([-cell for cell in chain])

def faceOrientation(boundaryLoop,face,FE,EV,cf):
    theBoundary = set(AA(ABS)(boundaryLoop))
    if theBoundary.intersection(FE[face])==set() and theBoundary.difference(FE[face])!=set(): ##BOH!!
        coboundaryFaces = [f for f in cf if set(FE[f]).intersection(theBoundary)!=set()]
        face = coboundaryFaces[0]            
    faceLoop,_ = triangulation.boundaryCycles(FE[face],EV)
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
    outputOp = larBoundary(*larModel)
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
        #VIEW(STRUCT( MKPOLS((V,[EV[h] for f in cf for h in FE[f]])) )) #add EV!
    return cf,coe

""" Main procedure of arrangement partitioning """
import inters,triangulation

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
    
    V,CV,FV,EV,CF,CE,COE = facesFromComponents((Z,FZ,EZ),FE,EF_angle)
    return V,CV,FV,EV,CF,CE,COE,FE

""" First steps of the Merge algorithm """    
def partition(W,FW,EW):
    quadArray = [[W[v] for v in face] for face in FW]
    parts = boxBuckets3d(containmentBoxes(quadArray))
    Z,FZ,EZ = spacePartition(W,FW,EW, parts)
    Z,FZ,EZ = inters.larSimplify((Z,FZ,EZ),radius=0.001)
    model = Z,FZ,EZ
    return model

def UBoundary2(EV,FV):
   UB2 = (coo2Csr(brc2Coo(EV)) * coo2Csr(brc2Coo(FV)).T).tocoo()
   data,i,j = TRANS([(1,i,j) for (d,i,j) in zip(UB2.data,UB2.row,UB2.col) if d==2])
   return coo_matrix((data,(i,j))).tocsc()

def SBoundary2(EV,FV):
   SB_2 = UBoundary2(EV,FV)
   for f in range(len(FV)):
      ind = defaultdict(list)
      chain_1 = list(SB_2[:,f].tocoo().row)
      chain_0 = zeros((2,len(chain_1)),dtype=int)
      for h,e in enumerate(chain_1):
         v1 = EV[e][0]
         v2 = EV[e][1]
         chain_0[0,h] = v1
         chain_0[1,h] = v2
         ind[v1] += [h]
         ind[v2] += [h]
      # tracking of ordered 0-chains
      k = 0
      while True:
         v1,v2 = chain_0[:,k]
         k = set(ind[v2]).difference([k]).pop()
         if k == 0: break
         if chain_0[0,k] == v2:
            v1,v2 = chain_0[:,k]
         else:
            chain_0[0,k],chain_0[1,k] = chain_0[1,k],chain_0[0,k]
      # TODO: 0-chains with multiple cycles ...
      # sign computation
      sign = []
      for h,e in enumerate(chain_1):
         v1,v2 = chain_0[:,h]
         if v1 < v2:
            sign += [1]
         else:
            sign += [-1]
      # update SB_2
      for h,e in enumerate(chain_1):
         SB_2[e,f] = sign[h]
   return SB_2

""" Next steps of the Merge algorithm """  
import boundary
  
def chain2coords(chain,n):
   data,i,j = TRANS([(code,cell,0) for (cell,code) in chain])
   coordVect = coo_matrix((data,(i,j)),(n,1))
   return coordVect
   
def coords2chainDict(coords):
   return dict(boundary.coords2chain(coords))

def next(cyclicPerm):
   def next1(pivot):
      ind = cyclicPerm.index(pivot)
      nextIndex = (ind + 1) % len(cyclicPerm)
      return cyclicPerm[nextIndex]
   return next1

def prev(cyclicPerm):
   def prev1(pivot):
      ind = cyclicPerm.index(pivot)
      nextIndex = (ind - 1) % len(cyclicPerm)
      return cyclicPerm[nextIndex]
   return prev1

def Choose(marks):
   try: return marks.index(1)
   except ValueError: return marks.index(0)

def SBoundary3(W,EW,FW):
   SB_2 = SBoundary2(EW,FW)
   D,I,J = [],[],[]
   m,n = SB_2.shape
   marks = [0 for k in range(n)]
   store = [0 for k in range(n)]
   # permutation subgroups of edges
   FE = [list(SB_2[:,f].tocoo().row) for f in range(SB_2.shape[1])]
   EF_angle, ET,TW,FT = faceSlopeOrdering((W,FW,EW),FE)
   cellCount = -1

   while sum(marks) < 2*n:
      # choose f
      f = Choose(marks)
      
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

""" Final steps of the Merge algorithm """  
def Boundary3(W,EW,FW):
   SB_3 = SBoundary3(W,EW,FW)
   vertDict = dict([(tuple(v),k) for k,v in enumerate(W)])
   vdict = sorted(vertDict)
   first, last = vertDict[vdict[0]], vertDict[vdict[-1]]
   WF = invertRelation(FW)
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

def MKSOLID(W,FW,EW):
   SB_2 = SBoundary2(EW,FW)
   FE = [list(SB_2[:,f].tocoo().row) for f in range(SB_2.shape[1])]
   triangleSet = boundaryTriangulation(W,FW,EW,FE)
   TW,FT = triangleIndices(triangleSet,W)
   B_3 = Boundary3(W,EW,FW)
   CF = [list(B_3[:,c].tocoo().row) for c in range(B_3.shape[1])]
   CT = [CAT([FT[f] for f in cell]) for cell in CF] 
   cells = AA(COMP([SOLIDIFY,STRUCT,MKPOLS]))(DISTL([ W, 
      [[TW[t] for t in cell] for cell in CT] ]))
   VIEW(EXPLODE(1.5,1.5,1.5)(cells))

