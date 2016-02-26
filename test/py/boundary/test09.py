""" Adjacency relations examples """
from larlib import *

sys.path.insert(0, 'test/py/boundary/')
from test07 import *

""" kfaces-to-kfaces relation """
eeOp = larEdges2Edges(EV,VV)
EE = [eeOp([k]) for k in range(len(EV))]

ffOp = larFaces2Faces(FV,EV)
FF = [ffOp([k]) for k in range(len(FV))]

ccOp = larCells2Cells(CV,FV,EV)
CC = [ccOp([k]) for k in range(len(CV))]


print "\nCC =",CC
print "\nFF =",FF
print "\nEE =",EE,"\n"

V,BF,BE = larBoundary3(V,CV,FV,EV)([1,0])
VIEW(STRUCT(MKTRIANGLES((V,[FV[h] for h in FF[-1]],EV),color=True)))
VIEW(STRUCT(MKPOLS((V,[EV[h] for h in EE[-1]]))+[COLOR(RED)(MKPOLS((V,[EV[-1]]))[0])]))
