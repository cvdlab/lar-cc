""" Visualize a list of colored HPCs for the faces in FV """
from larlib import *
    
filename = "test/svg/inters/graph3.svg"
lines = svg2lines(filename)
V,EV = lines2lar(lines)
VIEW(STRUCT(MKPOLS((V,EV))))
V,FV,EV,polygons = larFromLines(lines)

VV = AA(LIST)(range(len(V)))
hpc = STRUCT(MKPOLS((V,EV)))        
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],hpc,0.05))

VIEW(STRUCT(MKPOLYGONS(V,polygons)))
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLYGONS(V,polygons)))
VIEW(SKEL_1(STRUCT(MKPOLYGONS(V,polygons))))

viewLarComplexChain((V,FV,EV))([1])  # TODO: solve bug
