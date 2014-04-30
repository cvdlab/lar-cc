""" D LAR model input and handling """
from pyplasm import *
from scipy import *
import os,sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from lar2psm import *
from simplexn import *
from larcc import *
from largrid import *
from mapper import *
from boolean import vertexSieve

from architectural import *
from boolean import *

V = [[3,-3],
[9,-3],[0,0],[3,0],[9,0],[15,0],
[3,3],[6,3],[9,3],[15,3],[21,3], 
[0,9],[6,9],[15,9],[18,9],[0,13],
[6,13],[9,13],[15,13],[18,10],[21,10], 
[18,13],[6,16],[9,16],[9,17],[15,17],
[18,17],[-3,24],[6,24],[15,24],[-3,13]]
FV = [
[22,23,24,25,29,28], [15,16,22,28,27,30], [18,21,26,25], 
[13,14,19,21,18], [16,17,23,22], [11,12,16,15],
[9,10,20,19,14,13], [2,3,6,7,12,11], [0,1,4,8,7,6,3],
[4,5,9,13,18,17,16,12,7,8],[17,18,25,24,23]]
dwelling = [V,FV]

V,FV = larApply(t(3,0))(dwelling)
print "\n V,FV =",V,FV
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS(dwelling)))
dwelling = Struct([ t(3,0), dwelling ])
V1 = [[0,0],[3,0],[3,4.5],[0,4.5],[3,9],[0,9],[3,13],[-3,13],[-3,0],[0,-3]]
FV1 = [[0,1,2,3],[3,2,4,5],[0,3,5,4,6,7,8,9]]
landing = V1,FV1
stair = spiralStair(width=0.2,R=1.5,r=0.25,riser=0.1,pitch=4.,nturns=1.75,steps=36)
stair = larApply(r(0,0,PI/4))(stair)
# stair = larDisk([1.5,1.75*PI])()
# stair = larApply(r(PI/4))(stair)
stair = larApply(t(0,-3))(larCircle(2.5)())
plan = Struct([stair,landing,dwelling,s(-1,1),dwelling])
assembly2D = evalStruct(plan)
VIEW(EXPLODE(1.2,1.2,1)(CAT(AA(MKPOLS)(assembly2D))))
assembly1D = TRANS(CONS([S1,COMP([AA(face2edge),S2])])(TRANS(assembly2D)))
VIEW(EXPLODE(1.2,1.2,1)(CAT(AA(MKPOLS)(assembly1D))))
