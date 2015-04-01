""" Example of nested structures with translation and rotations """
import sys; sys.path.insert(0, 'lib/py/')
from largrid import *
from larstruct import *
square = larCuboids([1,1])
#square = Model(square)
table = larApply( t(-.5,-.5) )(square)
chair = Struct([ t(.75, 0), s(.35,.35), table ])
struct = Struct( [t(2,1)] + [table] + 4*[r(PI/2), chair])
struct = Struct(10*[struct,t(0,2.5)])
struct = Struct(10*[struct,t(3,0)])
scene = evalStruct(struct)
VIEW(SKEL_1(STRUCT(CAT(AA(MKPOLS)(scene)))))
