# prerequisite: 

from larlib import *
V,[VV,EV,FV,CV] = larCuboids([2,2,2],True)
cubeGrid = Struct([(V,FV,EV)],"cubeGrid")
cubeGrids = Struct(2*[cubeGrid,t(.25,.25,.25),r(0,0,PI/12)])
V,FV,EV = struct2lar(cubeGrids)
VV = AA(LIST)(range(len(V)))
frame = STRUCT(MKPOLS((V,EV)))
#VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV],frame,0.4))
quadArray = [[V[v] for v in face] for face in FV]
parts = boxBuckets3d(containmentBoxes(quadArray))

# faceFragmentation algorithm

def submodel(V,FV,EV):
	FE = crossRelation(FV,EV)
	def submodel0(f,F):
		fE = list(set(FE[f] + CAT([FE[g] for g in F])))
		fF = [f]+F
		return fF,fE
	return submodel0

def meetZero( sW, (w1,w2) ):
	testValue = sW[w1][2] * sW[w2][2]
	if testValue > 10**-4: 
		return False
	else: return True

def segmentIntersection(p1,p2):
	(x1,y1,z1),(x2,y2,z2) = p1,p2
	if abs(z1-z2) != 0.0:
		alpha = z1/(z1-z2)
		x = x1+alpha*(x2-x1)
		y = y1+alpha*(y2-y1)
		return x,y,0.0
	else: return None

def spacePartition(V,FV,EV, parts):
	submodel0 = submodel(V,FV,EV)
	out = []
	""" input: face index f; candidate incident faces F[f]; """
	for f,F in enumerate(parts):
		""" Selection of the LAR submodel S(f) := (V,FV,EV)(f) restricted to [f]+F[f] """	
		fF,fE = submodel0(f,F)
		subModel = Struct([(V,[FV[g] for g in fF],[EV[h] for h in fE])])
		sV,sFV,sEV = struct2lar(subModel)
		""" Computation of submanifold map M moving f to z=0 """
		pivotFace = [V[v] for v in FV[f]]
		M = submanifoldMapping(pivotFace)  # submanifold transformation
		""" Transformation of S(f) by M, giving S = (sW,sEW) := M(S(f)) """
		sW,sFW,sEW = larApply(M)((sV,sFV,sEV))
		""" filtering of EW edges traversing z=0, giving EZ edges and incident faces FZEZ """
		sFE = crossRelation(sFW,sEW)	
		edges = list(set([ e for k,face in enumerate(sFW)  for e in sFE[k] 
					if meetZero(sW, sEW[e]) ]))
		edgesPerFace = [ [e for e in sFE[k] if meetZero(sW, sEW[e])] 
					for k,face in enumerate(sFW) ]
		edges = list(set(CAT(edgesPerFace)))
		WW = AA(LIST)(range(len(sW)))
		wires = [sEW[e] for e in edges]
		wireFrame = STRUCT(MKPOLS((sW,wires)))
		""" for each face in FZEZ, computation of the aligned set of points p(z=0) """
		points = OrderedDict()
		lines = [[sW[w1],sW[w2]] for w1,w2 in wires]
		for k,(p,q) in enumerate(lines): 
			point = segmentIntersection(p,q)
			if point != None: points[edges[k]] = point
		pointsPerFace = [set(face).intersection(points.keys()) for face in edgesPerFace]
		lines = [[points[e][:2] for e in face] for face in pointsPerFace]
		lines = [line for line in lines if line!=[]]
		vpoints = [[(vcode(point),k) for k,point in enumerate(line)] for line in lines]
		lines = [AA(eval)(dict(line).keys()) for line in vpoints]
		### sorting of every aligned set FX, where X is the parametric coordinate along the intersection line
		### Check that every set FX has even cardinality
		### Construction of the planar set EX of lines including  EV and alternating lines from every set FX 
		""" (V2D,FV,EV) := larFromLines(EX) """
		z,fz,ez = larFromLines(lines)
		w,fw,ew = larApply(M.I)(([v+[0.0] for v in z],fz,ez))
		out += [Struct([(w,fw,ew)])]
		#VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((w,ew))))
	return struct2lar(Struct(out))

Z,FZ,EZ = spacePartition(V,FV,EV, parts)
VIEW(EXPLODE(1.1,1.1,1.1)(MKPOLS((Z,FZ,EZ))))