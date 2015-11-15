""" Generating the LAR of a set of non-intersecting cycles """
from larlib import *

sys.path.insert(0, 'test/py/inters/')
from test15 import *

cells = cellsFromCycles(testArray)
CV = AA(COMP([list,set,CAT]))(EVs)
EVdict = dict(zip(EV,range(len(EV))))
FE = [[EVdict[edge] for edge in cycle] for cycle in EVs] 
edges = [CAT([FE[cycle] for cycle in cell]) for cell in cells]
FVs = [[CV[cycle] for cycle in cell] for cell in cells]
FV = AA(CAT)(FVs)

def allBinarySubsetsOfLenght(n):
   out = [list(('{0:0'+str(n)+'b}').format(k)) for k in range(1,2**n)]
   return AA(AA(int))(out)

n = len(cells)
chains = allBinarySubsetsOfLenght(n)
"""
chains = [[1,0,0]]
chains += [[0,1,0]]
chains += [[0,0,1]]
chains += [[1,1,0]]
chains += [[0,1,1]]
chains += [[1,0,1]]
chains += [[1,1,1]]
"""
cycles = STRUCT(MKPOLS((V,EV)))
for chain in chains:
    chainBoundary = COLOR(RED)(STRUCT(MKPOLS((V,[EV[e] 
                        for e in chain2BoundaryChain(FV,EV)(chain)]))))
    VIEW(STRUCT([cycles, chainBoundary]))
