""" Splitting of othogonal lines """
import sys
sys.path.insert(0, 'lib/py/')
from inters import *

lines = [
[[0,0],[6,0]],
[[0,4],[6,4]],
[[0,0],[0,4]],
[[3,0],[3,4]],
[[6,0],[6,4]],
[[3,2],[6,2]]
]
VIEW(EXPLODE(1.2,1.2,1)(AA(POLYLINE)(lines)))

V,EV = lines2lar(lines)
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,EV))))
