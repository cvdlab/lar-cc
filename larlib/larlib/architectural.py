""" architectural module """
""" Initial import of modules """
from larlib import *

def face2edge(FV):
    """ From faces to list of edges """
    edges = AA(sorted)(CAT([TRANS([face, face[1:]+[face[0]]]) for face in FV]))
    return AA(eval)(set(AA(str)(edges)))
def lar2polylines (model):
    """ From LAR model to list of polylines """
    V,FV = model
    return [[V[v] for v in cell]+[V[cell[0]]] for cell in FV]
def lar2lines (model):
    """ From LAR model to list of lines """
    V,EV = model
    return [[V[v] for v in cell] for cell in EV]


def larCells(fun):
    def larCells0(assembly):
        return TRANS(CONS([S1,COMP([AA(fun),S2])])(TRANS(assembly)))
    return larCells0

def larVerts(fun):
    def larCells0(assembly):
        return TRANS(CONS([ COMP([AA(fun),S1]),S2 ])(TRANS(assembly)))
    return larCells0

def larBinOps(op):
    def larCells0(assembly):
        def larCells1(arg):
            return AA(op)(DISTR([assembly,arg]))
        return larCells1
    return larCells0


def larQuote1D(pattern):
    return larExtrude1( VOID, pattern )

def larQuote0D(pattern):
    V,CV = larQuote1D(pattern)
    return V,[[k] for k in range(len(V))] 

def bUnit_to_eEiP(FV,EV):
    """ Subdivide the 1-cells. 
    Return external envelope and interior partitions """
    eE = lar2boundaryEdges(FV,EV)
    iP = lar2InteriorEdges(FV,EV)
    return eE,iP

def lar2boundaryEdges(FV,EV):
    """ Boundary cells computation """
    return boundaryCells(FV,EV)

def lar2InteriorEdges(FV,EV):
    """ Boundary cells computation """
    boundarychain1 = boundaryCells(FV,EV)
    totalChain1 = range(len(EV))
    interiorCells = set(totalChain1).difference(boundarychain1)
    return interiorCells

def movePoint2point(twoModels):
    """ Move (P -> Q) operator """
    def movePoint2point0(pointP):
        def movePoint2point1(pointQ):
            [V,CV], [W,CW] = twoModels
            mat = t( *DIFF([pointP,pointQ]) )
            [W,CW] = larApply(mat)([W,CW])
            print "\n W =",W
            print "\n CW =",CW
            n = len(V)
            return [ V+W, CV+[[w+n for w in REVERSE(cell)] for cell in CW] ] 
        return movePoint2point1    
    return movePoint2point0


def spiralStair(width=0.2,R=1.,r=0.5,riser=0.1,pitch=2.,nturns=2.,steps=18):
    V,CV = larSolidHelicoid(width,R,r,pitch,nturns,steps)()
    W = CAT([[V[k],V[k+1],V[k+2],V[k+3]]+
        [SUM([V[k+1],[0,0,-riser]]),SUM([V[k+3],[0,0,-riser]])]
        for k,v in enumerate(V[:-4]) if k%4==0])
    for k,w in enumerate(W[:-12]):
        if k%6==0: W[k+1][2] = W[k+10][2]; W[k+3][2] = W[k+11][2]
    nsteps = len(W)/12
    CW =[SUM([[0,1,2,3,6,8,10,11],[6*k]*8]) for k in range(nsteps)]
    return W,CW

if __name__=="__main__":
    VIEW(STRUCT(MKPOLS(spiralStair())))
    VIEW(SKEL_1(STRUCT(MKPOLS(spiralStair()))))
    VIEW(STRUCT(MKPOLS(spiralStair(0.1))))

""" Solidify horizontal polygons in 3D """
def solidify(pol):    
    min=MIN([1])(pol)[0]
    max=MAX([1])(pol)[0]
    z = MIN([3])(pol)[0]
    pol = PROJECT(1)(pol)
    siz=max-min
    far_point=max+siz*100 
    def InftyProject(pol):
        verts,cells,pols=UKPOL(pol)
        verts=[[far_point] + v[1:] for v in verts]
        return MKPOL([verts,cells,pols])  
    ret=SPLITCELLS(pol)
    ret=[JOIN([pol,InftyProject(pol)]) for pol in ret]
    return T(3)(z)(XOR(ret))

def horizontalClosures(pattern):
    vertical1D = larQuote1D(pattern)
    vertical0D = larQuote0D(pattern)
    def horizontalClosures0(assembly2D):
        out = []
        for flat2D in assembly2D:
            V,FV = flat2D
            EV = face2edge(FV)
            dwellH3D = larModelProduct([flat2D,vertical0D])
            dwellV3D = larModelProduct([(V,EV),vertical1D])
            poly3D = AA(POLYLINE)(lar2polylines(dwellH3D))
            floors3D = AA(solidify)(poly3D)
            out += floors3D+MKPOLS(dwellV3D)
        return out
    return horizontalClosures0

""" Placing a 3D object (wall) with possible solid subtraction (door) """

def place(obj):

    def dist(p1,p2):
        return SQRT(SQR(p1[0]-p2[0])+SQR(p1[1]-p2[1]))

    depth,length,height = SIZE([1,2,3])(obj)
    p = array(MIN([1,2,3])(obj))
    obj = T([1,2,3])(list(-1*p))(obj)
    obj = S(2)(1./length)(obj)
    def place0(obj2=None):
        def place01(line):
            x,y,z = VECTDIFF([line[1],line[0]])
            angle = -math.atan2(x,y)
            outObj = S(2)(dist(line[1],line[0]))(obj)
            if isinstance(obj2,pyplasm.xgepy.Hpc):
                outObj = DIFFERENCE([outObj,obj2])
            outObj = R([1,2])(angle)(outObj)
            outObj = T([1,2,3])(line[0])(outObj)
            return outObj
        return place01
    return place0

