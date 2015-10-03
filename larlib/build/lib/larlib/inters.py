""" Module for pipelined intersection of geometric objects """
from larlib import *
from scipy import mat
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
def centroid(boxes,coord):
    delta,n = 0,len(boxes)
    ncoords = len(boxes[0])/2
    a = coord%ncoords
    b = a+ncoords
    for box in boxes:
        delta += (box[a] + box[b])/2
    return delta/n


""" XOR of FAN of ordered points """ 
def FAN(points): 
   pairs = zip(points[1:-2],points[2:-1])
   triangles = [MKPOL([[points[0],p1,p2],[[1,2,3]],None]) for p1,p2 in pairs]
   return XOR(triangles)
 
if __name__=="__main__":
    pol = [[0.476,0.332],[0.461,0.359],[0.491,0.375],[0.512,0.375],[0.514,0.375],
    [0.527,0.375],[0.543,0.34],[0.551,0.321],[0.605,0.314],[0.602,0.307],[0.589,
    0.279],[0.565,0.244],[0.559,0.235],[0.553,0.227],[0.527,0.239],[0.476,0.332]]

    VIEW(EXPLODE(1.2,1.2,1)(FAN(pol)))


""" Generation of random lines """
def randomLines(numberOfLines=200,scaling=0.3):
    randomLineArray = [redge(scaling) for k in range(numberOfLines)]
    [xs,ys] = TRANS(CAT(randomLineArray))[:2]
    xmin, ymin = min(xs), min(ys)
    v = array([-xmin,-ymin])
    randomLineArray = [[list(v1[:2]+v), list(v2[:2]+v)] for v1,v2 in randomLineArray]
    return randomLineArray

""" Containment boxes """
def containment2DBoxes(randomLineArray):
    boxes = [eval(vcode([min(x1,x2),min(y1,y2),max(x1,x2),max(y1,y2)]))
            for ((x1,y1),(x2,y2)) in randomLineArray]
    return boxes

""" Splitting the input above and below a threshold """
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


""" Iterative splitting of box buckets """
def splitting(bucket,below,above, finalBuckets,splittingStack):
    if (len(below)<4 and len(above)<4) or len(set(bucket).difference(below))<7 \
        or len(set(bucket).difference(above))<7: 
        finalBuckets.append(below)
        finalBuckets.append(above)
    else: 
        splittingStack.append(below)
        splittingStack.append(above)

def geomPartitionate(boxes,buckets):
    geomInters = [set() for h in range(len(boxes))]
    for bucket in buckets:
        for k in bucket:
            geomInters[k] = geomInters[k].union(bucket)
    for h,inters in enumerate(geomInters):
        geomInters[h] = geomInters[h].difference([h])
    return AA(list)(geomInters)

def boxBuckets(boxes):
    bucket = range(len(boxes))
    splittingStack = [bucket]
    finalBuckets = []
    while splittingStack != []:
        bucket = splittingStack.pop()
        below,above = splitOnThreshold(boxes,bucket,1)
        below1,above1 = splitOnThreshold(boxes,above,2)
        below2,above2 = splitOnThreshold(boxes,below,2)                      
        splitting(above,below1,above1, finalBuckets,splittingStack)
        splitting(below,below2,above2, finalBuckets,splittingStack)      
        finalBuckets = list(set(AA(tuple)(finalBuckets)))
    parts = geomPartitionate(boxes,finalBuckets)
    return AA(sorted)(parts)
    #return finalBuckets

""" Intersection of two line segments """
def segmentIntersect(boxes,lineArray,pointStorage):
    def segmentIntersect0(h):
        p1,p2 = lineArray[h]
        line1 = '['+ vcode(p1) +','+ vcode(p2) +']'
        (x1,y1),(x2,y2) = p1,p2
        B1,B2,B3,B4 = boxes[h]
        def segmentIntersect1(k):
            p3,p4 = lineArray[k]
            line2 = '['+ vcode(p3) +','+ vcode(p4) +']'
            (x3,y3),(x4,y4) = p3,p4
            b1,b2,b3,b4 = boxes[k]
            if not (b3<B1 or B3<b1 or b4<B2 or B4<b2):
            #if True:
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
                    if -0.0<=alpha<=1 and -0.0<=beta<=1:
                        pointStorage[line1] += [alpha]
                        pointStorage[line2] += [beta]
                        return list(array(p1)+alpha*(array(p2)-array(p1)))
            return None
        return segmentIntersect1
    return segmentIntersect0

""" Brute force bucket intersection """
def lineBucketIntersect(boxes,lineArray, h,bucket, pointStorage):
    intersect0 = segmentIntersect(boxes,lineArray,pointStorage)
    intersectionPoints = []
    intersect1 = intersect0(h)
    for line in bucket:
        point = intersect1(line)
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

    boxes = containment2DBoxes(lineArray)
    buckets = boxBuckets(boxes)
    intersectionPoints = set()
    for h,bucket in enumerate(buckets):
        pointBucket = lineBucketIntersect(boxes,lineArray, h,bucket, pointStorage)
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
    
    model = (V,EV)
    return larSimplify(model)

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

def facesFromComps(model):
    V,EV = model
    # Remove zero edges
    EV = list(set([ tuple(sorted([v1,v2])) for v1,v2 in EV if v1!=v2 ]))
    FV = []
    VE_angle = edgeSlopeOrdering((V,EV))
    csrEV = ordered_csrVE(VE_angle).T
    visited = zeros((len(EV),2))
    edge,v = firstSearch(visited)
    vertex = EV[edge][v]
    fv = []
    while True:
        if (edge,v) == (-1,-1):
            break #return [face for face in FV if face != None]
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
        FV = [face for face in FV if face != None]
    return V,FV,EV

""" SVG input parsing and transformation """
from larlib import *
import re # regular expression

def svg2lines(filename,containmentBox=[],rect2lines=True):
    stringLines = [line.strip() for line in open(filename)]   
    
    # SVG <line> primitives
    lines = [string.strip() for string in stringLines if re.match("<line ",string)!=None]   
    outLines = ""   
    for line in lines:
        searchObj = re.search( r'(<line )(.+)(" x1=")(.+)(" y1=")(.+)(" x2=")(.+)(" y2=")(.+)("/>)', line)
        if searchObj:
            outLines += "[["+searchObj.group(4)+","+searchObj.group(6)+"], ["+searchObj.group(8) +","+ searchObj.group(10) +"]],"
    if lines != []:
        lines = list(eval(outLines))
              
    # SVG <rect> primitives
    rects = [string.strip() for string in stringLines if re.match("<rect ",string)!=None]   
    outRects,searchObj = "",False 
    for rect in rects:
        searchObj = re.search( r'(<rect x=")(.+?)(" y=")(.+?)(" )(.*?)( width=")(.+?)(" height=")(.+?)("/>)', rect)
        if searchObj:
            outRects += "[["+searchObj.group(2)+","+searchObj.group(4)+"], ["+searchObj.group(8)+","+searchObj.group(10)+"]],"
    
    if rects != []:
        rects = list(eval(outRects))
        if rect2lines:
            lines += CAT([[[[x,y],[x+w,y]],[[x+w,y],[x+w,y+h]],[[x+w,y+h],[x,y+h]],[[x,y+h],[x,y]]] for [x,y],[w,h] in rects])
        else: 
            lines += [[[x,y],[x+w,y+h]] for [x,y],[w,h] in rects]
    for line in lines: print line
    
    """ SVG input normalization transformation """
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
    
    containmentBox = box
    
    return lines


""" Transformation of an array of lines in a 2D LAR complex """
def larFromLines(lines):
    V,EV = lines2lar(lines)
    V,EVs = biconnectedComponent((V,EV))
    EV = list(set(AA(tuple)(AA(sorted)(max(EVs, key=len)))))  ## NB
    V,EV = larRemoveVertices(V,EV)
    V,FV,EV = facesFromComps((V,EV))
    areas = integr.surfIntegration((V,FV,EV))
    boundaryArea = max(areas)
    interiorFaces = [FV[f] for f,area in enumerate(areas) if area!=boundaryArea and len(areas)>2]
    boundaryFace = FV[areas.index(boundaryArea)]
    return V,interiorFaces+[boundaryFace],EV

""" Pruning away clusters of close vertices """
from scipy.spatial import cKDTree

def pruneVertices(pts,radius=0.001):
   tree = cKDTree(pts)
   a = cKDTree.sparse_distance_matrix(tree,tree,radius)
   print a.keys()
   close = list(set(AA(tuple)(AA(sorted)(a.keys()))))
   import networkx as nx
   G=nx.Graph()
   G.add_nodes_from(range(len(pts)))
   G.add_edges_from(close)
   clusters, k, h = [], 0, 0
   
   subgraphs = list(nx.connected_component_subgraphs(G))
   V = [None for subgraph in subgraphs]
   vmap = [None for k in xrange(len(pts))]
   for k,subgraph in enumerate(subgraphs):
      group = subgraph.nodes()
      if len(group)>1: 
         V[k] = CCOMB([pts[v] for v in group])
         for v in group: vmap[v] = k
         clusters += [group]
      else: 
         oldNode = group[0]
         V[k] = pts[oldNode]
         vmap[oldNode] = k
   return V,close,clusters,vmap

""" Return a simplified LAR model """
def larSimplify(model,radius=0.001):
   if len(model)==2: V,CV = model 
   elif len(model)==3: V,CV,FV = model 
   else: print "ERROR: model input"
   
   W,close,clusters,vmap = pruneVertices(V,radius)
   celldim = DIM(MKPOL([V,[[v+1 for v in CV[0]]],None]))
   newCV = [list(set([vmap[v] for v in cell])) for cell in CV]
   CV = list(set([tuple(cell) for cell in newCV if len(cell) >= celldim+1]))
   CV = sorted(CV,key=len) # to get the boundary cell as last one (in most cases)

   if len(model)==3:
      celldim = DIM(MKPOL([V,[[v+1 for v in FV[0]]],None]))
      newFV = [list(set([vmap[v] for v in facet])) for facet in FV]
      FV = [facet for facet in newFV if len(facet) >= celldim]
      return W,CV,FV
   else: return W,CV

