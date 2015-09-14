
""" Testing of point-in-polygon classification algorithm """
from larlib import *

sys.path.insert(0, 'test/py/inters/')
from test10 import *


""" Half-line crossing test """
def crossingTest(new,old,count,status):
    if status == 0:
        status = new
        count += 0.5
    else:
        if status == old: count += 0.5
        else: count -= 0.5
        status = 0
    return count,status

""" Tile codes computation """
def setTile(box):
    tiles = [[9,1,5],[8,0,4],[10,2,6]]
    b1,b2,b3,b4 = box
    def tileCode(point):
        x,y = point
        code = 0
        if y>b1: code=code|1
        if y<b2: code=code|2
        if x>b3: code=code|4
        if x<b4: code=code|8
        return code 
    return tileCode

""" Point in polygon classification """
def pointInPolygonClassification(p,pol):
    x,y = p
    V,EV = pol
    xmin,xmax,ymin,ymax = x,x,y,y
    tilecode = setTile([ymax,ymin,xmax,xmin])
    count,status = 0,0
    for k,edge in enumerate(EV):
        p1,p2 = V[edge[0]],V[edge[1]]
        (x1,y1),(x2,y2) = p1,p2
        c1,c2 = tilecode(p1),tilecode(p2)
        k,c_edge, c_un, c_int = k,c1^c2, c1|c2, c1&c2
        #print "k,c_edge, c_un, c_int =",k,c_edge, c_un, c_int
        
        if c_edge == 0 and c_un == 0: return "p_on"
        elif c_edge == 12 and c_un == c_edge: return "p_on"
        elif c_edge == 3:
            if c_int == 0: return "p_on"
            elif c_int == 4: count += 1
        elif c_edge == 15:
            x_int = ((y-y2)*(x1-x2)/(y1-y2))+x2 
            if x_int > x: count += 1
            elif x_int == x: return "p_on"
        elif c_edge == 13 and ((c1==4) or (c2==4)):
                count,status = crossingTest(1,2,count,status)
        elif c_edge == 14 and (c1==4) or (c2==4):
                count,status = crossingTest(2,1,count,status)
        elif c_edge == 7: count += 1
        elif c_edge == 11: count = count
        elif c_edge == 1:
            if c_int == 0: return "p_on"
            elif c_int == 4: count,status = crossingTest(1,2,count,status)
        elif c_edge == 2:
            if c_int == 0: return "p_on"
            elif c_int == 4: count,status = crossingTest(2,1,count,status)
        elif c_edge == 4 and c_un == c_edge: return "p_on"
        elif c_edge == 8 and c_un == c_edge: return "p_on"
        elif c_edge == 5:
            if (c1==0) or (c2==0): return "p_on"
            else: count,status = crossingTest(1,2,count,status)
        elif c_edge == 6:
            if (c1==0) or (c2==0): return "p_on"
            else: count,status = crossingTest(2,1,count,status)
        elif c_edge == 9 and ((c1==0) or (c2==0)): return "p_on"
        elif c_edge == 10 and ((c1==0) or (c2==0)): return "p_on"
        #print "count,p1,p2 =",count,p1,p2        
    if (round(count)%2)==1: return "p_in"
    else: return "p_out"



result = []
for k in range(10000):
    queryPoint = [random.random(),random.random()]
    inOut = pointInPolygonClassification(queryPoint,[V,EV])
    #print k,queryPoint,inOut
    if inOut=="p_in": result += [MK(queryPoint)]
    elif inOut=="p_out": result += [COLOR(RED)(MK(queryPoint))]

VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(STRUCT([
    POLYLINE([[queryPoint[0],0],[queryPoint[0],1]]),
    POLYLINE([[0,queryPoint[1]],[1,queryPoint[1]]]),
    larModelNumbering(1,1,1)(V,[VV,EV,FV[:-1]],submodel,0.4)
    ]+result))
