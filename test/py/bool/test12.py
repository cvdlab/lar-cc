
""" Generation of the edge permutation associated to the 1-boundary of a 2-chain """
import sys;sys.path.insert(0, 'lib/py/')
from bool2 import *
sys.path.insert(0, 'test/py/larcc/')
from test11 import *

C2 = csr_matrix((len(FV),1))
for i in [21,16,23,22, 2,3,4, 9,28,5]: C2[i,0] = 1
BD = boundary(FV,EV)
C1 = BD * C2
C_1 = [i for i in range(len(EV)) if ABS(C1[i,0]) == 1 ]
C_2 = [i for i in range(len(FV)) if C2[i,0] == 1 ]

VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[EV[k] for k in C_1] + [FV[k] for k in C_2]))))

sign,next = cycles2permutation(boundaryCicles(C_1, EV))
