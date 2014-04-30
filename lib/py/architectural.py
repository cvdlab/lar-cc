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

def face2edge(FV):
   """ From faces to list of edges """
   edges = AA(sorted)(CAT([TRANS([face, face[1:]+[face[0]]]) for face in FV]))
   return AA(eval)(set(AA(str)(edges)))
def lar2polylines (model):
   """ From LAR model to list of polylines """
   V,FV = model
   return [[V[v] for v in cell]+[V[cell[0]]] for cell in FV]

def bUnit_to_eEiP(FV,EV):
   """ Subdivide the 1-cells. 
   Return external envelope and interior partitions """
   eE = lar2boundaryEdges(FV,EV)
   iP = lar2InteriorEdges(FV,EV)
   return eE,iP

def lar2boundaryEdges(FV,EV):
   """ Boundary cells computation """
   return boundaryCells(FV,EV)

def lar2InteriorEdges(FV,EV):
   """ Boundary cells computation """
   boundarychain1 = boundaryCells(FV,EV)
   totalChain1 = range(len(EV))
   interiorCells = set(totalChain1).difference(boundarychain1)
   return interiorCells

def movePoint2point(twoModels):
   """ Move (P -> Q) operator """
   def movePoint2point0(pointP):
      def movePoint2point1(pointQ):
         [V,CV], [W,CW] = twoModels
         mat = t( *DIFF([pointP,pointQ]) )
         [W,CW] = larApply(mat)([W,CW])
         print "\n W =",W
         print "\n CW =",CW
         n = len(V)
         return [ V+W, CV+[[w+n for w in REVERSE(cell)] for cell in CW] ] 
      return movePoint2point1 
   return movePoint2point0


def spiralStair(R=1.,r=0.5,riser=0.1,pitch=2.,nturns=2.,steps=18):
   V,CV = larSolidHelicoid(0.2,R,r,pitch,nturns,steps)()
   W = CAT([[V[k],V[k+1],V[k+2],V[k+3]]+
      [SUM([V[k+1],[0,0,-riser]]),SUM([V[k+3],[0,0,-riser]])]
      for k,v in enumerate(V[:-4]) if k%4==0])
   for k,w in enumerate(W[:-12]):
      if k%6==0: W[k+1][2] = W[k+10][2]; W[k+3][2] = W[k+11][2]
   nsteps = len(W)/12
   CW =[SUM([[0,1,2,3,6,8,10,11],[6*k]*8]) for k in range(nsteps)]
   return W,CW
   
VIEW(STRUCT(MKPOLS(spiralStair())))

