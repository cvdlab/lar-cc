import sys; sys.path.insert(0, 'lib/py/')
from largrid import *

mod_1 = larSplit(1)(4), larGrid(4)(1)
squares = larModelProduct([mod_1,mod_1])
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(squares)))
cubes = larModelProduct([squares,mod_1])
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(cubes)))
