""" Remove double instances of cells (and the unused vertices) """
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from larcc import *
from mapper import evalStruct


""" Remove the double instances of cells """
cellDict = defaultdict(list)
for k,cell in enumerate(FW):
    cellDict[tuple(cell)] += [k]
FW = [list(key) for key in cellDict.keys() if len(cellDict[key])==1]

VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((W,FW))))

""" Remove the unused vertices """
print "len(W) =",len(W)
V,FV = larRemoveVertices(W,FW)
print "len(V) =",len(V)

