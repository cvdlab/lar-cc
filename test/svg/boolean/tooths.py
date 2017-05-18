from larlib import *

filename = "test/svg/boolean/tooths.lines"
lines = array(lines2lines(filename))

v1,v2 = lines[0]
tx,ty = v1
vx,vy = array(v2) - v1
alpha = math.atan2(vy,vx)
trasl = t(-tx,-ty,0)
rot = r(0,0,alpha)

v1s = ((lines[:,0] * mat([[1,0,0,0],[0,1,0,0]]) + [0.,0.,0.,1.]) * trasl * rot).tolist()
v2s = ((lines[:,1] * mat([[1,0,0,0],[0,1,0,0]]) + [0.,0.,0.,1.]) * trasl * rot).tolist()
theLines = zip([v[:2] for v in v1s], [v[:2] for v in v2s])
pols = STRUCT(AA(POLYLINE)(theLines))
VIEW(pols)


#	TEST 1
#	======
def test1(theLines):
	V,FV,EV,polygons = larFromLines(theLines,True)
	VV = AA(LIST)(range(len(V)))
	"""
	VIEW(EXPLODE(1.2,1.2,1)(AA(SKEL_1)(MKPOLS((V,EV)))))
	polylines = AA(POLYLINE)([CAT([[V[v] for v in cycle+[cycle[0]]] for cycle in pol]) for pol in polygons])
	VIEW(EXPLODE(1.2,1.2,1)(polylines))
	colors = [CYAN,MAGENTA,WHITE,RED,YELLOW,GRAY,GREEN,ORANGE,BLUE,PURPLE,BROWN,BLACK]
	components = [COLOR(colors[k%12])(face) for k,face in enumerate(MKFACES((V,FV,EV)))]
	VIEW(STRUCT(components))
	"""

	W = (array(V) - V[6]).tolist()
	print "\nW =",W
	submodel = STRUCT(MKPOLS((W,EV)))        
	VIEW(larModelNumbering(1,1,1)(W,[VV,EV,FV],submodel,.2))
	Z = (array(W) * (mat([[1,0,0,0],[0,1,0,0]]) + [0.,0.,0.,1.]) * r(PI/2,0,0) * r(0,0,PI/6)).tolist()
	Z = [z[:3] for z in Z]
	VIEW(STRUCT(MKPOLS((Z,EV))))
	print "\nW[6] =",W[6]
	print "W[11] =",W[11]
	return Z

#	TEST 2
#	======
def test2(theLines):
	V,FV,EV,polygons = larFromLines(theLines,False)
	VV = AA(LIST)(range(len(V)))
	"""
	VIEW(EXPLODE(1.2,1.2,1)(AA(SKEL_1)(MKPOLS((V,EV)))))
	polylines = AA(POLYLINE)([CAT([[V[v] for v in cycle+[cycle[0]]] for cycle in pol]) for pol in polygons])
	VIEW(EXPLODE(1.2,1.2,1)(polylines))
	colors = [CYAN,MAGENTA,WHITE,RED,YELLOW,GRAY,GREEN,ORANGE,BLUE,PURPLE,BROWN,BLACK]
	components = [COLOR(colors[k%12])(face) for k,face in enumerate(MKFACES((V,FV,EV)))]
	VIEW(STRUCT(components))
	"""

	W = (array(V) - V[16]).tolist()
	submodel = STRUCT(MKPOLS((W,EV)))        
	VIEW(larModelNumbering(1,1,1)(W,[VV,EV,FV],submodel,300))
	Z = (array(W) * (mat([[1,0,0,0],[0,1,0,0]]) + [0.,0.,0.,1.]) * r(PI/2,0,0) * r(0,0,PI/6)).tolist()
	Z = [z[:3] for z in Z]
	VIEW(STRUCT(MKPOLS((Z,EV))))
	print "\nW[16] =",W[16]
	print "W[13] =",W[13]
	return Z
