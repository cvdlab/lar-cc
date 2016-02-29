""" Incidence relations examples """
from larlib import *

sys.path.insert(0, 'test/py/boundary/')
from test08 import *

""" mfaces-to-nfaces relations """
fcOp = larCells2Faces(CV,FV,EV)
CF = [fcOp([k]) for k in range(len(CV))]
FC = invertRelation(CF)

ecOp = larCells2Edges(CV,FV,EV)
CE = [ecOp([k]) for k in range(len(CV))]
EC = invertRelation(CE)
    
efOp = larFaces2Edges(FV,EV)
FE = [efOp([k]) for k in range(len(FV))]
EF = invertRelation(FE)

print "\nFC =",FC
print "\nEC =",EC
print "\nEF =",EF,"\n"


