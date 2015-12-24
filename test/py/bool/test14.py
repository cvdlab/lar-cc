
""" Testing of point-in-polygon classification  on a random polygon """
from larlib import *

filename = "test/svg/inters/tile.svg"
lines = svg2lines(filename)
V,EV = lines2lar(lines)



result = []
classify = pointInPolygonClassification((V,EV))
for k in range(10000):
    queryPoint = [random.random(),random.random()]
    inOut = classify(queryPoint)
    #print k,queryPoint,inOut
    if inOut=="p_in": result += [MK(queryPoint)]
    elif inOut=="p_out": result += [COLOR(RED)(MK(queryPoint))]

VIEW(STRUCT(MKPOLS((V,EV))+result))
