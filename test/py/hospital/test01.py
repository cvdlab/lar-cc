import sys
sys.path.insert(0, 'lib/py/')
from hospital import *

models = evalStruct(copyStruct(groundFloor))
lars = AA(AS(SEL)([1,3]))(models)
VIEW(STRUCT(CAT(AA(MKPOLS)(lars))))

Vs,EVs = TRANS(AA(AS(SEL)([1,3]))(models))
Vs = AA(metric)(Vs)
lars = TRANS([Vs,EVs])
VIEW(STRUCT(CAT(AA(MKPOLS)(lars))))
