""" Example of oriented edge drawing """
import sys;sys.path.insert(0, 'lib/py/')
from larcc import *
sys.path.insert(0, 'test/py/larcc/')
from test11 import *

C2 = csr_matrix((len(FV),1))
for i in [1,2, 12,13,14,15, 22,23, 29,30,31]: C2[i,0] = 1
BD = boundary(FV,EV)
C1 = BD * C2
C_1 = [i for i in range(len(EV)) if ABS(C1[i,0]) == 1 ]
C_2 = [i for i in range(len(FV)) if C2[i,0] == 1 ]

VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[EV[k] for k in C_1] + [FV[k] for k in C_2]))))
