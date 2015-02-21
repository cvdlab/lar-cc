
""" Generate an array of 3D lines and their boxes """
from pyplasm import *
""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *

randomTriaArray = randomTriangles(8000,0.1)
VIEW(STRUCT(AA(POLYLINE)(randomTriaArray)))

boxes = containmentBoxes(randomTriaArray)
hexas = AA(box2exa)(boxes)
glass = MATERIAL([1,0,0,0.1,  0,1,0,0.1,  0,0,1,0.1, 0,0,0,0.1, 100])
cyan = COLOR(CYAN)(STRUCT(AA(POLYLINE)(randomTriaArray)))
yellow = STRUCT(AA(glass)(AA(MKPOL)(hexas)))
VIEW(STRUCT([cyan,yellow]))

buckets = boxBuckets(boxes)
colors = [CYAN, MAGENTA, WHITE, RED, YELLOW, GRAY, GREEN, ORANGE, BLACK, BLUE, PURPLE, BROWN]
sets = [COLOR(colors[k%12])(STRUCT(AA(POLYLINE)([randomTriaArray[h] 
            for h in bucket]))) for k,bucket in enumerate(buckets)]
VIEW(STRUCT(sets))
