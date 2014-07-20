""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')
from bool import *

V1 = [[0,0],[10,0],[10,10],[0,10]]
FV1 = [range(4)]
EV1 = [[0,1],[1,2],[2,3],[3,0]]
VV1 = AA(LIST)(range(len(V1)))

V2 = [[2.5,2.5],[7.5,2.5],[7.5,7.5],[2.5,7.5]]
FV2 = [range(4)]
EV2 = [[0,1],[1,2],[2,3],[3,0]]
VV2 = AA(LIST)(range(len(V2)))
""" Bulk of Boolean task computation """
""" Computation of edges an input visualisation """
model1 = V1,FV1
model2 = V2,FV2
basis1 = [VV1,EV1,FV1]
basis2 = [VV2,EV2,FV2]
submodel12 = STRUCT(MKPOLS((V1,EV1))+MKPOLS((V2,basis2[1])))
VIEW(larModelNumbering(V1,basis1,submodel12,4))
VIEW(larModelNumbering(V2,basis2,submodel12,4))


V,CV,chain1,chain2 = booleanChains((V1,basis1), (V2,basis2))

VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[cell for cell,c in zip(CV,chain1) if c==1] ))))
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[cell for cell,c in zip(CV,chain2) if c==1] ))))
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[cell for cell,c1,c2 in zip(CV,chain1,chain2) if c1+c2==2] ))))
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[cell for cell,c1,c2 in zip(CV,chain1,chain2) if c1+c2==1] ))))
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[cell for cell,c1,c2 in zip(CV,chain1,chain2) if c1+c2>=1] ))))


