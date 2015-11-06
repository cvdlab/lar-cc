
""" Testing of point-in-polygon classification algorithm """
from larlib import *

filename = "test/svg/inters/tile.svg"
lines = svg2lines(filename)
V,EV = lines2lar(lines)


""" Half-line crossing test """
def crossingTest(new,old,count,status):
    if status == 0:
        status = new
        count += 0.5
    else:
        if status == old: count += 0.5
        else: count -= 0.5
        status = 0

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
def pointInPolygonClassification(pol):

    V,EV = pol
    # edge orientation
    FV = [sorted(set(CAT(EV)))]
    orientedCycles = boundaryPolylines(Struct([(V,FV,EV)]))
    EV = []
    for cycle in orientedCycles:
        EV += zip(cycle[:-1],cycle[1:])

    def pointInPolygonClassification0(p):
        x,y = p
        xmin,xmax,ymin,ymax = x,x,y,y
        tilecode = setTile([ymax,ymin,xmax,xmin])
        count,status = 0,0
    
        for k,edge in enumerate(EV):
            p1,p2 = edge[0],edge[1]
            (x1,y1),(x2,y2) = p1,p2
            c1,c2 = tilecode(p1),tilecode(p2)
            c_edge, c_un, c_int = c1^c2, c1|c2, c1&c2
            
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
                    crossingTest(1,2,status,count)
            elif c_edge == 14 and (c1==4) or (c2==4):
                    crossingTest(2,1,status,count)
            elif c_edge == 7: count += 1
            elif c_edge == 11: count = count
            elif c_edge == 1:
                if c_int == 0: return "p_on"
                elif c_int == 4: crossingTest(1,2,status,count)
            elif c_edge == 2:
                if c_int == 0: return "p_on"
                elif c_int == 4: crossingTest(2,1,status,count)
            elif c_edge == 4 and c_un == c_edge: return "p_on"
            elif c_edge == 8 and c_un == c_edge: return "p_on"
            elif c_edge == 5:
                if (c1==0) or (c2==0): return "p_on"
                else: crossingTest(1,2,status,count)
            elif c_edge == 6:
                if (c1==0) or (c2==0): return "p_on"
                else: crossingTest(2,1,status,count)
            elif c_edge == 9 and ((c1==0) or (c2==0)): return "p_on"
            elif c_edge == 10 and ((c1==0) or (c2==0)): return "p_on"
        if ((round(count)%2)==1): return "p_in"
        else: return "p_out"
    return pointInPolygonClassification0



result = []
classify = pointInPolygonClassification((V,EV))
for k in range(10000):
    queryPoint = [random.random(),random.random()]
    inOut = classify(queryPoint)
    #print k,queryPoint,inOut
    if inOut=="p_in": result += [MK(queryPoint)]
    elif inOut=="p_out": result += [COLOR(RED)(MK(queryPoint))]

VIEW(STRUCT(MKPOLS((V,EV))+result))
