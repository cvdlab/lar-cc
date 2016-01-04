""" 2D polygon triangulation """
from larlib import *
#from support import PolygonTessellator,vertex

filename = "test/svg/bool/interior.svg"
lines = svg2lines(filename)    
V,FV,EV,polygons = larFromLines(lines)
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,FV[:-1]+EV)) + AA(MK)(V)))

VIEW(EXPLODE(1.2,1.2,1)(MKTRIANGLES((V,FV,EV)) ))
