""" Extraction of boundary vertices of a cuboidal complex """
import sys; sys.path.insert(0, 'lib/py/')
from largrid import *

shape = (10,10,10)
model = larCuboids(shape)
V,cells = model
exterior = cuboidalComplexBoundaryVertices(model)
VIEW(STRUCT(MKPOLS((V,AA(LIST)(exterior)))))
V,facets = larFacets((V,cells+[exterior]))
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(larFacets((V,cells+[exterior])))))
FV = improperFacetsCovering(facets,cells)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV))))
