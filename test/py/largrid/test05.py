""" Extraction of boundary vertices of a cuboidal complex """
import sys; sys.path.insert(0, 'lib/py/')
from largrid import *

shape = (50,50)
model = larCuboids(shape)
V,cells = model
exterior = cuboidalComplexBoundaryVertices(model)
VIEW(STRUCT(MKPOLS((V,AA(LIST)(exterior)))))
VIEW(STRUCT(MKPOLS(larFacets((V,cells),dim=2))))
V,facets = larFacets((V,cells+[exterior]),dim=2)
EV = improperFacetsCovering(facets,cells,2)
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((V,EV))))
