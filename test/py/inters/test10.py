""" Biconnected components from orthogonal LAR model """
import sys
sys.path.insert(0, 'lib/py/')
from inters import *
from bool1 import larRemoveVertices
from hospital import surfIntegration
from iot3d import polyline2lar

""" SVG input parsing and transformation """
from larcc import *
import re # regular expression

filename = "test/py/inters/test1.svg"
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
VIEW(STRUCT(AA(POLYLINE)(lines)))

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
VIEW(STRUCT(AA(POLYLINE)(lines)))


V,EV = lines2lar(lines)




print "\nV =",V
print "\nEV =",EV
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,EV))))
model = V,EV
VV = vertices2vertices(model)
leaves = [k for k,vv in enumerate(VV) if len(vv)==1]
EV_ = [[v1,v2]  for v1,v2 in EV if set(leaves).intersection([v1,v2]) == set()]
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,EV_))))

EV = list(set(AA(tuple)(sorted(AA(sorted)(EV_))))) 
V,EV = larRemoveVertices(V,EV)
VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV],submodel,0.10))

model = V,EV
FV = facesFromComponents((V,EV))
areas = surfIntegration((V,FV,EV))
boundaryArea = max(areas)
faces = [FV[f] for f,area in enumerate(areas) if area!=boundaryArea]
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,faces+EV)) + AA(MK)(V)))

V,FV,EV = polyline2lar([[ V[v] for v in FV[areas.index(boundaryArea)] ]])
VIEW(STRUCT(MKPOLS((V,EV))))
