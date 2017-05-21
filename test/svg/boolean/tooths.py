from larlib import *
	
#	TEST 1
#	======
def test1(lines):
	V,FV,EV,polygons = larFromLines(lines,True)
	V = (array(V) - V[14]).tolist()
	VV = AA(LIST)(range(len(V)))
	submodel = STRUCT(MKPOLS((V,EV)))        
	VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,.2)) 
	Z = (array(V) * (mat([[1,0,0,0],[0,1,0,0]]) + [0.,0.,0.,1.]) * r(PI/2,0,0).T * r(0,0,PI/6).T ).tolist()
	Z = [z[:3] for z in Z]
	VIEW(STRUCT(MKPOLS((Z,EV))))
	return Z,FV,EV

#	TEST 2
#	======
def test2(lines):
	V,FV,EV,polygons = larFromLines(lines,False)
	V = (array(V) - V[14]).tolist()
	VV = AA(LIST)(range(len(V)))
	submodel = STRUCT(MKPOLS((V,EV)))        
	VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,300))
	Z = (array(V) * (mat([[1,0,0,0],[0,1,0,0]]) + [0.,0.,0.,1.]) * r(PI/2,0,0).T * r(0,0,PI/6).T).tolist()
	Z = [z[:3] for z in Z]
	VIEW(STRUCT(MKPOLS((Z,EV))))
	return Z,FV,EV
	
	
filename = "test/svg/boolean/tooths.lines"
lines = lines2lines(filename)

v1,v2 = lines[0]
vs = array(CAT(lines)) - v1
vx,vy = vs[1]
alpha = -math.atan2(vy,vx)
M = mat([[cos(alpha),-sin(alpha)],[sin(alpha),cos(alpha)]])
W = (M * vs.T).T.tolist()
myLines = [[W[2*k],W[2*k+1]] for k in range(len(W)/2)]

V,FV,EV,polygons = larFromLines(myLines,False)
VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))        
VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,300))

W,FV,EV = test1(myLines)
W,FV,EV = test2(myLines)


