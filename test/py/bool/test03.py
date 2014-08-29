""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')
from bool import *

V1 = [[0,0],[10,0],[10,10],[0,10]]
FV1 = [range(4)]
EV1 = [[0,1],[1,2],[2,3],[0,3]]
VV1 = AA(LIST)(range(len(V1)))

V2 = [[2.5,2.5],[12.5,2.5],[12.5,12.5],[2.5,12.5]]
FV2 = [range(4)]
EV2 = [[0,1],[1,2],[2,3],[0,3]]
VV2 = AA(LIST)(range(len(V2)))
""" Bulk of Boolean task computation """
""" Computation of edges an input visualisation """
model1 = V1,FV1
model2 = V2,FV2
basis1 = [VV1,EV1,FV1]
basis2 = [VV2,EV2,FV2]
submodel12 = STRUCT(MKPOLS((V1,EV1))+MKPOLS((V2,EV2)))
VIEW(larModelNumbering(V1,basis1,submodel12,4))
VIEW(larModelNumbering(V2,basis2,submodel12,4))


V,CV,BC,CVbits,vertdict,dict_facets,splittingCovectors,n_bf1,n_bf2 = \
   booleanChains((V1,basis1), (V2,basis2))
   
""" Boundary-Coboundary operators in the SCDC basis """
dim = len(V[0])
FV = larConvexFacets (V,CV)
EV = larFacets((V,FV), dim)

VV = AA(LIST)(range(len(V)))
if dim == 3: bases = [VV,EV,FV,CV]
elif dim == 2: bases = [VV,FV,CV]
else: print "\nerror: not implemented\n"

coBoundaryMat = signedCellularBoundary(V,bases).T
""" Boundaries in SCDC coordinates """

def removeDuplicates(dictOfLists):
   values = AA(COMP([sorted,list,set,AA(abs)]))(dictOfLists.values())
   keys = dictOfLists.keys()
   return dict(zip(keys,values))

def splitBoundaries(CV,FV,n_bf1,n_bf2,splittingCovectors,coBoundaryMat):
   bucket1,bucket2 = defaultdict(list),defaultdict(list) 
   for k,cell in enumerate(CV):
      if splittingCovectors[k] != []:
         facets = list(coBoundaryMat[k].tocoo().col)
         signs = list(coBoundaryMat[k].tocoo().data)
         orientedFacets = AA(prod)(zip(facets,signs))
         for facet in orientedFacets:
            for face,covector,verts in splittingCovectors[k]:
               if all([ verySmall(INNERPROD([covector,V[v]+[1.0]])) 
                        for v in FV[abs(facet)] ]):
                  if face<n_bf1: bucket1[face] += [facet]
                  elif face<n_bf1+n_bf2: bucket2[face] += [facet]
                  else: print "error: separation of argument boundaries"
                  break    
   facets1 = removeDuplicates(bucket1)
   facets2 = removeDuplicates(bucket2)
   return facets1,facets2

""" Coboundary of boundaries """

FV = larConvexFacets (V,CV)
boundary1,boundary2 = splitBoundaries(CV,FV,n_bf1,n_bf2,splittingCovectors,
                              coBoundaryMat)
facets1 = [FV[facet] for face in boundary1 for facet in boundary1[face]]
facets2 = [FV[facet] for face in boundary2 for facet in boundary2[face]]
submodel1 = mkSignedEdges((V,facets1))
submodel2 = mkSignedEdges((V,facets2))

VIEW(STRUCT([submodel1,submodel2]))
VIEW(STRUCT([ COLOR(YELLOW)(EXPLODE(1.2,1.2,1)(MKPOLS((V,facets1)))), 
      COLOR(GREEN)(EXPLODE(1.2,1.2,1)(MKPOLS((V,facets2)))) ]))
      
boundaryMat = coBoundaryMat.T

def cellTagging(boundaryDict,boundaryMat,CV,FV,V,BC,CVbits,arg):
   dim = len(V[0])
   for face in boundaryDict:
      for facet in boundaryDict[face]:
         cofaces = list(boundaryMat[facet].tocoo().col)
         cosigns = list(boundaryMat[facet].tocoo().data)
         if len(cofaces) == 1: 
            CVbits[cofaces[0]][arg] = 1
         elif len(cofaces) == 2:
            v0 = list(set(CV[cofaces[0]]).difference(FV[facet]))[0]
            v1 = list(set(CV[cofaces[1]]).difference(FV[facet]))[0]
            # take d affinity independent vertices in face (TODO)
            simplex0 = BC[face][:dim] + [v0]
            simplex1 = BC[face][:dim] + [v1]
            sign0 = sign(det([V[v]+[1] for v in simplex0]))
            sign1 = sign(det([V[v]+[1] for v in simplex1]))
            if sign0 == 1: CVbits[cofaces[0]][arg] = 1
            elif sign0 == -1: CVbits[cofaces[0]][arg] = 0
            if sign1 == 1: CVbits[cofaces[1]][arg] = 1
            elif sign1 == -1: CVbits[cofaces[1]][arg] = 0
         else: 
            print "error: too many cofaces of boundary facets"
   return CVbits
   
CVbits = cellTagging(boundary1,boundaryMat,CV,FV,V,BC,CVbits,0)
CVbits = cellTagging(boundary2,boundaryMat,CV,FV,V,BC,CVbits,1)



VV = AA(LIST)(range(len(V)))
FV = larConvexFacets (V,CV)
submodel = STRUCT(MKPOLS((V,FV)))
VIEW(larModelNumbering(V,[VV,FV,CV],submodel,3))

chain1,chain2 = TRANS(CVbits)

if DEBUG:
   VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[cell for cell,c in zip(CV,chain1) if c==1] ))))
   VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[cell for cell,c in zip(CV,chain2) if c==1] ))))
   # VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[cell for cell,c1,c2 in zip(CV,chain1,chain2) if c1*c2==1] ))))
   VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[cell for cell,c1,c2 in zip(CV,chain1,chain2) if c1+c2==1] ))))
   VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,[cell for cell,c1,c2 in zip(CV,chain1,chain2) if c1+c2>=1] ))))

