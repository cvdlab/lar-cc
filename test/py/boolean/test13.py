""" Polygon triangulation and importing to LAR """

from larlib import *
from p2t import *

""" load points """
def load_points(file_name):
    infile = open(file_name, "r")
    points = []
    while infile:
        line = infile.readline()
        s = line.split()
        if len(s) == 0:
            break
        points.append([float(s[0]), float(s[1])])
    return points

""" input polyline visualization """
points = load_points("test/data/nazca_monkey.dat")
VIEW(POLYLINE(points))

""" CDT triangulation with poly2tri """
polyline = [Point(p[0],p[1]) for p in points]  
cdt = CDT(polyline)
triangles = cdt.triangulate()

""" conversion of triangulation to LAR """
trias = [ [[t.a.x,t.a.y],[t.b.x,t.b.y],[t.c.x,t.c.y]] for t in triangles ]
vdict = defaultdict(list)
for k,point in enumerate(CAT(trias)): 
    vdict[vcode(4)(point)] += [k]
vdict = OrderedDict(zip(vdict.keys(),range(len(vdict.keys()))))
FV = [(vdict[vcode(4)(a)], vdict[vcode(4)(b)], vdict[vcode(4)(c)]) for [a,b,c] in trias] 
repeatedEdges = CAT([[[v1,v2],[v2,v3],[v3,v1]] for [v1,v2,v3] in FV])
EV = list(set(AA(tuple)(AA(sorted)(repeatedEdges))))
V = [eval(vect) for vect in vdict]

""" visualization of LAR model of triangulation """
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,FV))))
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,EV))))

""" reconstruction of boundary polyline for LAR model """
EW = boundaryCells(FV,EV)
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[EV[e] for e in EW]))))
model = (V,FV,[EV[e] for e in EW])
struct = Struct([model])

""" visualization of LAR generated boundary polyline """
poly = boundaryModel2polylines(structBoundaryModel(struct))
VIEW(POLYLINE(poly[0][:-1]))

