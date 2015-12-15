""" Testing containments between non intersecting cycles """
from larlib import *

filename = "test/svg/inters/facade.svg"
lines = svg2lines(filename)
VIEW(STRUCT(AA(POLYLINE)(lines)))

V,EV = lines2lar(lines)
V,EVs = biconnectedComponent((V,EV))
# candidate face
FVs = AA(COMP([list,set,CAT]))(EVs)

testArray = latticeArray(V,EVs)

for k in range(len(testArray)):
   print k,testArray[k]
print "\ncells = ", cellsFromCycles(testArray),"\n"

VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FVs],submodel,0.15)) 
