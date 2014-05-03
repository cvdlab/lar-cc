""" 3D mock-up of apartment block """
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
plan = Struct([landing,dwelling,s(-1,1),dwelling])
assembly2D = evalStruct(plan)
assembly1D = larCells(face2edge)(assembly2D)
VIEW(EXPLODE(1.2,1.2,1)(CAT(AA(MKPOLS)(assembly1D))))

stair = spiralStair(width=0.2,R=3,r=0.25,riser=0.1,pitch=4.4,nturns=1.75,steps=36)
stair = larApply(r(0,0,3*PI/4))(stair)
stair = larApply(t(0,-3,0))(stair)
stairColumn = larApply(t(0,-3,0))(larRod(0.25,4.2)())
mod_1 = larQuote1D( 6*[0.2,-3.8] )
assembly3D = larBinOps(larModelProduct)(assembly2D)(mod_1)
VIEW(EXPLODE(1.2,1.2,1)(CAT(AA(MKPOLS)(assembly3D))))

horClosures = horizontalClosures([0.2,-3.8]*12 +[0.2])(assembly2D)
VIEW(STRUCT(horClosures))

wire = SKEL_1(INSR(PROD)(AA(QUOTE)([[6,9,9,9,9,6],[-3,10,11],[4]*12])))
VIEW(wire)
frame3D = T(1)(-24)(OFFSET([.2,.6,.2])(wire))
VIEW(frame3D)

assembly3D = evalStruct(Struct([stairColumn,stair,t(0,0,4)]*12))
VIEW(STRUCT(CAT(AA(MKPOLS)(assembly3D)) + horClosures + [frame3D]))
