""" Test example of LAR of a 2-complex with non-contractible and non-manifold cells"""
from larlib import *
sys.path.insert(0, 'test/py/triangulation/')
from test09 import *

model = V,FV,EV
faces = MKFACES(model)
colors = [CYAN,MAGENTA,WHITE,RED,YELLOW,GRAY,GREEN,ORANGE,BLACK,BLUE,PURPLE,BROWN]
components = [COLOR(colors[k%12])(face) for k,face in enumerate(faces)]
VIEW(STRUCT(components))
