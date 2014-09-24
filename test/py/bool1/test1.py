
import sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from bool1 import *

""" Definition of Boolean arguments """
V1 = [[3,0],[11,0],[13,10],[10,11],[8,11],[6,11],[4,11],[1,10],[4,3],[6,4],
      [8,4],[10,3]]
FV1 = [[0,1,8,9,10,11],[1,2,11],[3,10,11],[4,5,9,10],[6,8,9],[0,7,8],[2,3,
      11],[3,4,10],[5,6,9],[6,7,8]]
EV1 = [[0,1],[0,7],[0,8],[1,2],[1,11],[2,3],[2,11],[3,4],[3,10],[3,11],[4,
      5],[4,10],[5,6],[5,9],[6,7],[6,8],[6,9],[7,8],[8,9],[9,10],[10,11]]
VV1 = AA(LIST)(range(len(V1)))

V2 = [[0,3],[14,2],[14,5],[14,7],[14,11],[0,8],[3,7],[3,5]]
FV2 = [[0,5,6,7],[0,1,7],[4,5,6],[2,3,6,7],[1,2,7],[3,4,6]]
EV2 = [[0,1],[0,5],[0,7],[1,2],[1,7],[2,3],[2,7],[3,4],[3,6],[4,5],[4,6],
      [5,6],[6,7]]
VV2 = AA(LIST)(range(len(V2)))

arg1 = V1,(VV1,EV1,FV1)
arg2 = V2,(VV2,EV2,FV2)

def makeFacetDicts(FW,boundary1,boundary2):
   FWdict = dict()
   for k,facet in enumerate (FW): FWdict[str(facet)] = k
   for key,value in boundary1.items():
      value = [FWdict[str(facet)] for facet in value]
      boundary1[key] = value
   for key,value in boundary2.items():
      value = [FWdict[str(facet)] for facet in value]
      boundary2[key] = value
   return boundary1,boundary2,FWdict

def larBool(arg1,arg2):
   V1,basis1 = arg1
   V2,basis2 = arg2
   cells1 = basis1[-1]
   cells2 = basis2[-1]
   model1,model2 = (V1,cells1),(V2,cells2)
   
   """ First Boolean step """
   def larBool1():
      V, CV1,CV2, n1,n12,n2 = mergeVertices(model1,model2)
      VV = AA(LIST)(range(len(V)))
      V,CV,vertDict,n1,n12,n2,BC,nbc1,nbc2 = makeCDC(arg1,arg2)
      W,CW,VC,BCellCovering,cellCuts,boundary1,boundary2,BCW = makeSCDC(V,CV,BC,nbc1,nbc2)
      assert len(VC) == len(V) 
      assert len(BCellCovering) == len(BC)
      return W,CW,VC,BCellCovering,cellCuts,boundary1,boundary2,BCW 
   
   """ Second Boolean step """
   def larBool2(boundary1,boundary2):
      dim = len(W[0])
      WW = AA(LIST)(range(len(W)))
      FW = larConvexFacets (W,CW)
      if len(CW)==4: FW=[[0,1],[1,2],[0,2],[0,3],[2,3],[2,4],[2,5],
                     [3,4],[4,5]] #test5.py
      _,EW = larFacets((W,FW), dim=2)
      boundary1,boundary2,FWdict = makeFacetDicts(FW,boundary1,boundary2)
      if dim == 3: 
         _,EW = larFacets((W,FW), dim=2)
         bases = [WW,EW,FW,CW]
      elif dim == 2: bases = [WW,FW,CW]
      else: print "\nerror: not implemented\n"
      return W,CW,dim,bases,boundary1,boundary2,FW,BCW
   
   """ Third Boolean step """
   def larBool3():
      coBoundaryMat = signedCellularBoundary(W,bases).T
      boundaryMat = coBoundaryMat.T
      CWbits = [[-1,-1] for k in range(len(CW))]
      CWbits = cellTagging(boundary1,boundaryMat,CW,FW,W,BCW,CWbits,0)
      CWbits = cellTagging(boundary2,boundaryMat,CW,FW,W,BCW,CWbits,1)
      for cell in range(len(CW)):
         if CWbits[cell][0] == 1:
            CWbits = booleanChainTraverse(0,cell,W,CW,CWbits,1)      
         if CWbits[cell][0] == 0:
            CWbits = booleanChainTraverse(0,cell,W,CW,CWbits,0)
         if CWbits[cell][1] == 1:
            CWbits = booleanChainTraverse(1,cell,W,CW,CWbits,1)
         if CWbits[cell][1] == 0:
            CWbits = booleanChainTraverse(1,cell,W,CW,CWbits,0)
      chain1,chain2 = TRANS(CWbits)
      return W,CW,FW,boundaryMat,boundary1,boundary2,chain1,chain2   
   
   """ Fourth Boolean step """
   def larBool4(W):
      W,CX = gatherPolytopes(W,CW,FW,boundaryMat,boundary1,boundary2)
      FX = larConvexFacets (W,CX)      
      return W,CX,FX
   
      
   W,CW,VC,BCellCovering,cellCuts,boundary1,boundary2,BCW = larBool1()
   W,CW,dim,bases,boundary1,boundary2,FW,BCW = larBool2(boundary1,boundary2)
   W,CW,FW,boundaryMat,boundary1,boundary2,chain1,chain2 = larBool3()
   W,CX,FX = larBool4(W)

   def larBool0(op): 
      if op == "union":
         chain = [cell for cell,c1,c2 in zip(CW,chain1,chain2) if c1+c2>=1]
      elif op == "intersection":
         chain = [cell for cell,c1,c2 in zip(CW,chain1,chain2) if c1*c2==1]
      elif op == "xor":
         chain = [cell for cell,c1,c2 in zip(CW,chain1,chain2) if c1+c2==1]
      elif op == "difference":
         chain = [cell for cell,c1,c2 in zip(CW,chain1,chain2) if c1==1 and c2==0 ]
      else: print "Error: non implemented op"
      return W,CW,chain,CX,FX
   return larBool0
   
boolean = larBool(arg1,arg2)  

W,CW,chain,CX,FX = boolean("xor")
VIEW(EXPLODE(1.1,1.1,1)(MKPOLS((W,chain))))
W,CW,chain,CX,FX = boolean("union")
VIEW(EXPLODE(1.1,1.1,1)(MKPOLS((W,chain))))
W,CW,chain,CX,FX = boolean("intersection")
VIEW(EXPLODE(1.1,1.1,1)(MKPOLS((W,chain))))
W,CW,chain,CX,FX = boolean("difference")
VIEW(EXPLODE(1.1,1.1,1)(MKPOLS((W,chain))))

VIEW(EXPLODE(1.1,1.1,1)(MKPOLS((W,CX))))
VIEW(SKEL_1(EXPLODE(1.1,1.1,1)(MKPOLS((W,FX)))))

