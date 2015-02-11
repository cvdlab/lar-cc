""" Module for pipelined intersection of geometric objects """
from pyplasm import *
""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')
from larcc import *
DEBUG = True

""" Coding utilities """
""" Generation of a random point """
def rpoint():
    return eval( vcode([ random.random(), random.random() ]) )

""" Generation of a random line segment """
def redge(scaling):
    v1,v2 = array(rpoint()), array(rpoint())
    c = (v1+v2)/2
    pos = rpoint()
    v1 = (v1-c)*scaling + pos
    v2 = (v2-c)*scaling + pos
    return tuple(eval(vcode(v1))), tuple(eval(vcode(v2)))

""" Transformation of a 2D box into a closed polyline """    
def box2rect(box):
    x1,y1,x2,y2 = box
    verts = [[x1,y1],[x2,y1],[x2,y2],[x1,y2],[x1,y1]]
    return verts

""" Computation of the 1D centroid of a list of 2D boxes """    
def centroid(boxes,xy='x'):
    delta,n = 0,len(boxes)
    if xy=='x': a=0; b=2
    elif xy=='y': a=1; b=3
    for box in boxes:
        delta += (box[a] + box[b])/2
    return delta/n


""" Generation of random lines """
def randomLines(numberOfLines=200,scaling=0.3):
    randomLineArray = [redge(scaling) for k in range(numberOfLines)]
    [xs,ys] = TRANS(CAT(randomLineArray))
    xmin, ymin = min(xs), min(ys)
    v = array([-xmin,-ymin])
    randomLineArray = [[list(v1+v), list(v2+v)] for v1,v2 in randomLineArray]
    return randomLineArray

""" Containment boxes """
def containmentBoxes(randomLineArray):
    boxes = [eval(vcode([min(x1,x2),min(y1,y2),max(x1,x2),max(y1,y2)]))
            for ((x1,y1),(x2,y2)) in randomLineArray]
    return boxes

""" Splitting the input above and below a threshold """
def splitOnThreshold(boxes,subset,xy='x'):
    theBoxes = [boxes[k] for k in subset]
    threshold = centroid(theBoxes,xy)
    if xy=='x': a=0;b=2;
    elif xy=='y': a=1;b=3;
    below,above = [],[]
    for k in subset:
        if boxes[k][a] <= threshold: below += [k]
    for k in subset:
        if boxes[k][b] >= threshold: above += [k]
    return below,above


""" Iterative splitting of box buckets """
def boxBuckets(boxes):
    bucket = range(len(boxes))
    splittingStack = [bucket]
    finalBuckets = []
    while splittingStack != []:
        bucket = splittingStack.pop()
        below,above = splitOnThreshold(boxes,bucket,'x')
        below1,above1 = splitOnThreshold(boxes,above,'y')
        below2,above2 = splitOnThreshold(boxes,below,'y')
        
        if (len(below1)<4 and len(above1)<4) or len(set(bucket).difference(below1))<7 \
            or len(set(bucket).difference(above1))<7: 
            finalBuckets.append(below1)
            finalBuckets.append(above1)
        else: 
            splittingStack.append(below1)
            splittingStack.append(above1)
            
        if (len(below2)<4 and len(above2)<4) or len(set(bucket).difference(below2))<7 \
            or len(set(bucket).difference(above2))<7:  
            finalBuckets.append(below2)
            finalBuckets.append(above2)
        else: 
            splittingStack.append(below2)
            splittingStack.append(above2)
    return list(set(AA(tuple)(finalBuckets)))

""" Intersection of two line segments """
def segmentIntersect(pointStorage):
    def segmentIntersect0(segment1):
        p1,p2 = segment1
        line1 = '['+ vcode(p1) +','+ vcode(p2) +']'
        (x1,y1),(x2,y2) = p1,p2
        #B1,B2,B3,B4 = eval(vcode([min(x1,x2),min(y1,y2),max(x1,x2),max(y1,y2)]))
        def segmentIntersect1(segment2):
            p3,p4 = segment2
            line2 = '['+ vcode(p3) +','+ vcode(p4) +']'
            (x3,y3),(x4,y4) = p3,p4
            #b1,b2,b3,b4 = eval(vcode([min(x3,x4),min(y3,y4),max(x3,x4),max(y3,y4)]))
            #if ((B1<=b1<=B3) or (B1<=b3<=B3)) and ((B2<=b2<=B4) or (B2<=b4<=B4)):
            if True:
                m23 = mat([p2,p3])
                m14 = mat([p1,p4])
                m = m23 - m14
                v3 = mat([p3])
                v1 = mat([p1])
                v = v3-v1
                a=m[0,0]; b=m[0,1]; c=m[1,0]; d=m[1,1];
                det = a*d-b*c
                if det != 0:
                    m_inv = mat([[d,-b],[-c,a]])*(1./det)
                    alpha, beta = (v*m_inv).tolist()[0]
                    #alpha, beta = (v*m.I).tolist()[0]
                    if 0<=alpha<=1 and 0<=beta<=1:
                        pointStorage[line1] += [alpha]
                        pointStorage[line2] += [beta]
                        return list(array(p1)+alpha*(array(p2)-array(p1)))
            return None
        return segmentIntersect1
    return segmentIntersect0

""" Brute force bucket intersection """
def lineBucketIntersect(lines,pointStorage):
    intersect0 = segmentIntersect(pointStorage)
    intersectionPoints = []
    n = len(lines)
    for k,line in enumerate(lines):
        intersect1 = intersect0(line)
        for h in range(k+1,n):
            line1 = lines[h]
            point = intersect1(line1)
            if point != None: 
                intersectionPoints.append(eval(vcode(point)))
    return intersectionPoints

""" Accelerate intersection of lines """
def lineIntersection(lineArray):

    from collections import defaultdict
    pointStorage = defaultdict(list)
    for line in lineArray:
        p1,p2 = line
        key = '['+ vcode(p1) +','+ vcode(p2) +']'
        pointStorage[key] = []

    boxes = containmentBoxes(lineArray)
    buckets = boxBuckets(boxes)
    intersectionPoints = set()
    for bucket in buckets:
        lines = [lineArray[k] for k in bucket]
        pointBucket = lineBucketIntersect(lines,pointStorage)
        intersectionPoints = intersectionPoints.union(AA(tuple)(pointBucket))

    frags = AA(eval)(pointStorage.keys())
    params = AA(COMP([sorted,list,set,tuple,eval,vcode]))(pointStorage.values())
        
    return intersectionPoints,params,frags  ### GOOD: 1, WRONG: 2 !!!

""" Create the LAR of fragmented lines """
from scipy import spatial

def lines2lar(lineArray):
    _,params,frags = lineIntersection(lineArray)
    vertDict = dict()
    index,defaultValue,V,EV = -1,-1,[],[]
    
    for k,(p1,p2) in enumerate(frags):
        outline = [vcode(p1)]
        if params[k] != []:
            for alpha in params[k]:
                if alpha != 0.0 and alpha != 1.0:
                    p = list(array(p1)+alpha*(array(p2)-array(p1)))
                    outline += [vcode(p)]
        outline += [vcode(p2)]
    
        edge = []
        for key in outline:
            if vertDict.get(key,defaultValue) == defaultValue:
                index += 1
                vertDict[key] = index
                edge += [index]
                V += [eval(key)]
            else:
                edge += [vertDict[key]]
            EV.extend([[edge[k],edge[k+1]] for k,v in enumerate(edge[:-1])])
    
    # identification of close vertices
    closePairs = scipy.spatial.KDTree(V).query_pairs(10**(-PRECISION))
    if closePairs != []:
        EV_ = []
        for v1,v2 in EV:
            for v,w in closePairs:
                if v1 == w: v1 = v
                elif v2 == w: v2 = v
            EV_ += [[v1,v2]]
        EV = EV_
        print "\nclosePairs =",closePairs

    # Remove double edges
    EV = list(set(AA(tuple)(AA(sorted)(EV))))

    return V,EV

""" Biconnected components """
""" Adjacency lists of 1-complex vertices """
def vertices2vertices(model):
    V,EV = model
    csrEV = csrCreate(EV)
    csrVE = csrTranspose(csrEV)
    csrVV = matrixProduct(csrVE,csrEV)    
    cooVV = csrVV.tocoo()
    data,rows,cols = AA(list)([cooVV.data, cooVV.row, cooVV.col])
    triples = zip(data,rows,cols)
    VV = [[] for k in range(len(V))]
    for datum,row,col in triples:
        if row != col: VV[col] += [row]
    return AA(sorted)(VV)

""" Main procedure for biconnected components """
def biconnectedComponent(model):
    W,_ = model
    V = range(len(W))
    count = 0
    stack,out = [],[]
    visited = [None for v in V]
    parent = [None for v in V]
    d = [None for v in V]
    low = [None for v in V]
    for u in V: visited[u] = False
    for u in V: parent[u] = []
    VV = vertices2vertices(model)
    for u in V: 
        if not visited[u]: 
            DFV_visit( VV,out,count,visited,parent,d,low,stack, u )
    return W,[component for component in out if len(component) > 1]

""" Hopcroft-Tarjan algorithm """
def DFV_visit( VV,out,count,visited,parent,d,low,stack,u ):
    visited[u] = True
    count += 1
    d[u] = count
    low[u] = d[u]
    for v in VV[u]:
        if not visited[v]:
            stack += [(u,v)]
            parent[v] = u
            DFV_visit( VV,out,count,visited,parent,d,low,stack, v )
            if low[v] >= d[u]:
                out += [outputComp(stack,u,v)]
            low[u] = min( low[u], low[v] )
        else:
            if not (parent[u]==v) and (d[v] < d[u]):
                stack += [(u,v)]
                low[u] = min( low[u], d[v] )

""" Output of biconnected components """
def outputComp(stack,u,v):
    out = []
    while True:
        e = stack.pop()
        out += [list(e)]
        if e == (u,v): break
    return list(set(AA(tuple)(AA(sorted)(out))))


""" Circular ordering of edges around vertices """
def edgeSlopeOrdering(model):
    V,EV = model
    from bool1 import invertRelation
    VE,VE_angle = invertRelation(EV),[]
    for v,ve in enumerate(VE):
        ve_angle = []
        if ve != []:
            for edge in ve:
                v0,v1 = EV[edge]
                if v == v0:     x,y = list(array(V[v1]) - array(V[v0]))
                elif v == v1:    x,y = list(array(V[v0]) - array(V[v1]))
                angle = math.atan2(y,x)
                ve_angle += [180*angle/PI]
        pairs = sorted(zip(ve_angle,ve))
        #VE_angle += [TRANS(pairs)[1]]
        VE_angle += [[pair[1] for pair in pairs]]
    return VE_angle

""" Ordered incidence relationship of vertices and edges """
def ordered_csrVE(VE_angle):
    triples = []
    for v,ve in enumerate(VE_angle):
        n = len(ve)
        for k,edge in enumerate(ve):
            triples += [[v, ve[k], ve[ (k+1)%n ]]]
    csrVE = triples2mat(triples,shape="csr")
    return csrVE

""" Faces from biconnected components """

def firstSearch(visited):
    for edge,vertices in enumerate(visited):
        for v,vertex in enumerate(vertices):
            if visited[edge,v] == 0.0:
                visited[edge,v] = 1.0
                return edge,v
    return -1,-1

def facesFromComponents(model):
    V,EV = model
    FV = []
    VE_angle = edgeSlopeOrdering(model)
    csrEV = ordered_csrVE(VE_angle).T
    visited = zeros((len(EV),2))
    edge,v = firstSearch(visited)
    vertex = EV[edge][v]
    fv = []
    while True:
        if (edge,v) == (-1,-1):
            return [face for face in FV if face != None]
        elif (fv == []) or (fv[0] != vertex):
            
            fv += [vertex]
            nextEdge = csrEV[edge,vertex]
            v0,v1 = EV[nextEdge]
            
            try:
                vertex, = set([v0,v1]).difference([vertex])
            except ValueError:
                print 'ValueError: too many values to unpack'
                break
                
            if v0==vertex: pos=0
            elif v1==vertex: pos=1
                        
            if visited[nextEdge, pos] == 0:
                visited[nextEdge, pos] = 1
                edge = nextEdge                
        else:
            FV += [fv]
            fv = []
            edge,v = firstSearch(visited)
            vertex = EV[edge][v]
        #print "fv =",fv
        #print "edge,vertex =",edge,vertex
    return [face for face in FV if face != None]

""" SVG input parsing and transformation """
from larcc import *
import re # regular expression

def svg2lines(filename):

    lines = [line.strip() for line in open(filename) if re.match("<line ",line)!=None]   
    for line in lines: print line
        
    out = ""    
    for line in lines:
        #searchObj = re.search( r'([0-9]*\.[0-9]*)(.*?)([0-9]*\.[0-9]*)(.*?)([0-9]*\.[0-9]*)(.*?)([0-9]*\.[0-9]*)', line)
        searchObj = re.search( r'(<line )(.+)(" x1=")(.+)(" y1=")(.+)(" x2=")(.+)(" y2=")(.+)("/>)', line)
        if searchObj:
            #out += "[["+searchObj.group(1)+","+searchObj.group(3)+"], ["+searchObj.group(5)+","+searchObj.group(7)+"]],"
            out += "[["+searchObj.group(4)+","+searchObj.group(6)+"], ["+searchObj.group(8)+","+searchObj.group(10)+"]],"
    
    lines = list(eval(out))
    
    # window-viewport transformation
    xs,ys = TRANS(CAT(lines))
    box = [min(xs), min(ys), max(xs), max(ys)]
    
    # viewport aspect-ratio checking, setting a computed-viewport 'b'
    b = [None for k in range(4)]
    if (box[2]-box[0])/(box[3]-box[1]) > 1:  
        b[0]=0; b[2]=1; bm=(box[3]-box[1])/(box[2]-box[0]); b[1]=.5-bm/2; b[3]=.5+bm/2
    else: 
        b[1]=0; b[3]=1; bm=(box[2]-box[0])/(box[3]-box[1]); b[0]=.5-bm/2; b[2]=.5+bm/2
    
    # isomorphic 'box -> b' transform to standard unit square
    lines = [[[ 
    ((x1-box[0])*(b[2]-b[0]))/(box[2]-box[0]) , 
    ((y1-box[1])*(b[3]-b[1]))/(box[1]-box[3]) + 1], [
    ((x2-box[0])*(b[2]-b[0]))/(box[2]-box[0]), 
    ((y2-box[1])*(b[3]-b[1]))/(box[1]-box[3]) + 1]]  
          for [[x1,y1],[x2,y2]] in lines]
    
    # line vertices set to fixed resolution
    lines = eval("".join(['['+ vcode(p1) +','+ vcode(p2) +'], ' for p1,p2 in lines]))
    return lines

