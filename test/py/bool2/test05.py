""" Face (and incident faces) transformation """

""" Two unit cubes """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *
glass = MATERIAL([1,0,0,0.1,  0,1,0,0.1,  0,0,1,0.1, 0,0,0,0.1, 100])

V,[VV,EV,FV,CV] = larCuboids([1,1,1],True)
cube1 = Struct([(V,FV,EV)],"cube1")
twoCubes = Struct([cube1,t(.5,.5,.5),cube1])
V,FV,EV = struct2lar(twoCubes)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV))))

quadArray = [[V[v] for v in face] for face in FV]
boxes = containmentBoxes(quadArray)
hexas = AA(box2exa)(boxes)
parts = boxBuckets(boxes)


for k,bundledFaces in enumerate(parts):
    edges,faces = intersection(V,FV,EV)(k,bundledFaces)
    for face in faces:
        line = []
        for edge in face:
            (x1,y1,z1),(x2,y2,z2) = edge
            if not verySmall(z2-z1):
                x = (x2-x1)/(z2-z1) + x1
                y = (y2-y1)/(z2-z1) + y1
                p = [x,y,0]
                line += [eval(vcode(p))]
        edges += [line]
    print "k,edges =",k,edges
    
    hpcedges = AA(POLYLINE)(edges)
    VIEW(STRUCT(hpcedges))
    v,fv,ev = larFromLines([[point[:-1] for point in edge] for edge in edges])
    VIEW(EXPLODE(1.2,1.2,1)( MKPOLS((v,ev)) ))
    
    vv = AA(LIST)(range(len(v)))
    submodel = STRUCT(MKPOLS((v,ev)))
    VIEW(larModelNumbering(1,1,1)(v,[vv,ev,fv[:-1]],submodel,1))
    
    polylines = [[v[k] for k in face+[face[0]]] for face in fv[:-1]]
    VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((v,ev)) + AA(MK)(v) + AA(FAN)(polylines) ))



