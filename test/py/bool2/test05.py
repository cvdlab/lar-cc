
""" Preparation of test data for the Boolean algorithm """
from pyplasm import *
""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *
from mapper import *

mod1 = larSplit(1)(1), larGrid(1)(1)
mod_1 = struct2lar(Struct(5*[mod1,t(2)]))
squares = larModelProduct([mod_1,mod_1])
#VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(squares)))
cubes = larModelProduct([squares,mod_1])
obj1 = struct2lar(Struct([s(1./5,1./5,1./5), r(PI/12,PI/12,PI/12), t(-4.5,-4.5,-4.5), cubes]))

basis = larDisk(0.5)([18,5])
height = larSplit(1)(6), larGrid(6)(1)
obj2 = larModelProduct([basis,height])

VIEW(STRUCT(CAT(AA(MKPOLS)([obj1,obj2]))))
boxes = lar2boxes(obj1,1)+lar2boxes(obj2,2)
hexas = TRANS(AA(box2exa)(boxes))[0]
types = TRANS(AA(box2exa)(boxes))[1]
VIEW(EXPLODE(2,2,2)(AA(MKPOL)(hexas)))

buckets = boxBuckets(boxes)
colors = [CYAN, MAGENTA, WHITE, RED, YELLOW, GRAY, GREEN, ORANGE, BLACK, BLUE, PURPLE, BROWN]
target = []
for k,bucket in enumerate(buckets):
    if set([types[h][0] for h in bucket])=={1,2}:
        print "\nk,types =",k,[types[h] for h in bucket]
        VIEW(STRUCT(AA(view(obj1,obj2))([types[h] for h in bucket])))
        target += [COLOR(colors[k%12])( MKPOL(hexas[h]) ) ]
VIEW(EXPLODE(1.2,1.2,1.2)(target))
