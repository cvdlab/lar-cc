""" Generation of independent line buckets """
import sys
sys.path.insert(0, 'lib/py/')
from inters import *

lines = randomLines(200,0.3)
VIEW(STRUCT(AA(POLYLINE)(lines)))

boxes = containmentBoxes(lines)
buckets = boxBuckets(boxes)

colors = [CYAN, MAGENTA, WHITE, RED, YELLOW, GRAY, GREEN, ORANGE, BLACK, BLUE, PURPLE, BROWN]
sets = [COLOR(colors[k%12])(STRUCT(AA(POLYLINE)([lines[h] 
            for h in bucket]))) for k,bucket in enumerate(buckets)]

VIEW(STRUCT(sets))
