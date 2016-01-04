""" Generating the LAR of a set of non-intersecting cycles """
from larlib import *

sys.path.insert(0, 'test/py/triangulation/')
from test01 import *

cells = cellsFromCycles(latticeArray)
CV = AA(COMP([list,set,CAT]))(EVs)
EVdict = dict(zip(EV,range(len(EV))))
FE = [[EVdict[edge] for edge in cycle] for cycle in EVs] 
edges = [CAT([FE[cycle] for cycle in cell]) for cell in cells]
FVs = [[CV[cycle] for cycle in cell] for cell in cells]
FV = AA(CAT)(FVs)
cycles = STRUCT(MKPOLS((V,EV)))
csrBoundaryMat = boundary(FV,EV)

n = len(cells)
chains = allBinarySubsetsOfLenght(n)
for chain in chains:
    chainBoundary = COLOR(RED)(STRUCT(MKPOLS((V,[EV[e] 
                        for e in chain2BoundaryChain(csrBoundaryMat)(chain)]))))
    VIEW(STRUCT([cycles, chainBoundary]))
