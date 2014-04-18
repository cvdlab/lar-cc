
from pyplasm import *

def larSplit(dom):
    def larSplit1(n):
        # assert n > 0 and isinstance(n,int)
        item = float(dom)/n
        ints = range(n+1)
        items = [item]*(n+1)
        vertices = [[int*item] for (int,item) in zip(ints,items)]
        return vertices
    return larSplit1

assert larSplit(1)(3) == [[0.0], [0.3333333333333333], [0.6666666666666666], [1.0]]
assert larSplit(1)(1) == [[0.0], [1.0]]
assert larSplit(2*PI)(12) == [[0.0], [0.5235987755982988], [1.0471975511965976], 
[1.5707963267948966], [2.0943951023931953], [2.617993877991494], 
[3.141592653589793], [3.665191429188092], [4.1887902047863905], 
[4.71238898038469], [5.235987755982988], [5.759586531581287], 
[6.283185307179586]]
