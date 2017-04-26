""" Biconnected components from orthogonal LAR model """
from larlib import *
colors = [CYAN, MAGENTA, YELLOW, RED, GREEN, ORANGE, PURPLE, WHITE, BLACK, BLUE]

lines = randomLines(300,.4)
V,EV = lines2lar(lines)
model = V,EV
VIEW(STRUCT(AA(POLYLINE)(lines)))


boxes = containment2DBoxes(lines)
rects= AA(box2rect)(boxes)
cyan = COLOR(CYAN)(STRUCT(AA(POLYLINE)(lines)))
yellow = COLOR(YELLOW)(STRUCT(AA(POLYLINE)(rects)))
VIEW(STRUCT([cyan,yellow]))



V,EVs = biconnectedComponent(model)
HPCs = [STRUCT(MKPOLS((V,EV))) for EV in EVs]
sets = [COLOR(colors[k%10])(hpc) for k,hpc in enumerate(HPCs)]
VIEW(STRUCT(sets))

EV = CAT(EVs)
V,EV = larRemoveVertices(V,EV)
V,FV,EV = facesFromComps((V,EV))
areas = surfIntegration((V,FV,EV))
boundaryArea = max(areas)
FV = [FV[f] for f,area in enumerate(areas) if area!=boundaryArea]

polylines = [[V[v] for v in face+[face[0]]] for face in FV]
VIEW(EXPLODE(1.2,1.2,1)( AA(FAN)(polylines) + AA(COLOR(CYAN))(MKPOLS((V,EV))) + AA(COLOR(RED))(AA(MK)(V)) ))

colors = [CYAN, MAGENTA, WHITE, RED, YELLOW, GRAY, GREEN, ORANGE, BLACK, BLUE, PURPLE, BROWN]
sets = [COLOR(colors[k%12])(FAN(pol)) for k,pol in enumerate(polylines)]
VIEW(STRUCT(sets))


VIEW(EXPLODE(1.2,1.2,1)((AA(FAN)(polylines))))
VIEW(EXPLODE(1.2,1.2,1)((AA(POLYLINE)(polylines))))

VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
VIEW(larModelNumbering(1,1,1)(V,[VV,EV],submodel,0.1))
