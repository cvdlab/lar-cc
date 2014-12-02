""" Extraction of oriented boundary of a cuboidal 2-complex """
import sys; sys.path.insert(0, 'lib/py/')
from largrid import *
from random import random

shape = 20,20
V,cells = larCuboids(shape)
cellSpan = prod(shape)
fraction = 0.9
remove = [int(random()*cellSpan) for k in range(int(cellSpan*fraction)) ]
cells = [cells[k] for k in range(cellSpan) if not k in remove]
V,EV = larCuboidsFacets((V,cells))
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,EV))))

boundaryCells = signedBoundaryCells(V,cells,EV)
def swap(mylist): return [mylist[1]]+[mylist[0]]+mylist[2:]
orientedBoundary = [EV[-k] if k<0 else swap(EV[k]) for k in boundaryCells]
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,orientedBoundary))))

