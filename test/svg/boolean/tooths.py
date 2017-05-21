from larlib import *
	
#	TEST 1
#	======
def test1(lines):
	V,FV,EV,polygons = larFromLines(lines,True)
	VV = AA(LIST)(range(len(V)))
	submodel = STRUCT(MKPOLS((V,EV)))        
	VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,.2))
	Z = (array(V) * (mat([[1,0,0,0],[0,1,0,0]]) + [0.,0.,0.,1.]) * r(PI/2,0,0) * r(0,0,PI/6)).tolist()
	Z = [z[:3] for z in Z]
	VIEW(STRUCT(MKPOLS((Z,EV))))
	print "\nV[6] =",V[6]
	print "V[11] =",V[11]
	return Z

#	TEST 2
#	======
def test2(lines):
	V,FV,EV,polygons = larFromLines(lines,False)
	VV = AA(LIST)(range(len(V)))
	submodel = STRUCT(MKPOLS((V,EV)))        
	VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],submodel,.2))
	print "V =",V
	Z = (array(V) * (mat([[1,0,0,0],[0,1,0,0]]) + [0.,0.,0.,1.]) * r(PI/2,0,0) * r(0,0,PI/6)).tolist()
	Z = [z[:3] for z in Z]
	VIEW(STRUCT(MKPOLS((Z,EV))))
	print "\nW[16] =",V[16]
	print "W[13] =",V[13]
	return Z
	
	
filename = "test/svg/boolean/tooths.lines"
lines = lines2lines(filename)
myLines,N = normalize(lines)
mypol = STRUCT(AA(POLYLINE)(myLines))
VIEW(mypol)

v1,v2 = myLines[0]
vs = array(CAT(myLines)) - v1
vx,vy = vs[1]
alpha = -math.atan2(vy,vx)
M = mat([[cos(alpha),-sin(alpha)],[sin(alpha),cos(alpha)]])
W = (M * vs.T).T.tolist()
myLines = [[W[2*k],W[2*k+1]] for k in range(len(W)/2)]
mypol = STRUCT(AA(POLYLINE)(myLines))
VIEW(mypol)



test1(lines)
test2(lines)


