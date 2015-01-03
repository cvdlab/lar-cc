"""Module with functions needed to interface LAR with pyplasm"""
def importModule(moduleName):
    import sys; sys.path.insert(0, 'lib/py/')
    import moduleName
    

from lar2psm import *
def t(*args): 
    d = len(args)
    mat = scipy.identity(d+1)
    for k in range(d): 
        mat[k,d] = args[k]
    return mat.view(Mat)

def s(*args): 
    d = len(args)
    mat = scipy.identity(d+1)
    for k in range(d): 
        mat[k,k] = args[k]
    return mat.view(Mat)

def r(*args): 
    args = list(args)
    n = len(args)
    if n == 1: # rotation in 2D
        angle = args[0]; cos = COS(angle); sin = SIN(angle)
        mat = scipy.identity(3)
        mat[0,0] = cos;    mat[0,1] = -sin;
        mat[1,0] = sin;    mat[1,1] = cos;
    
    if n == 3: # rotation in 3D
        mat = scipy.identity(4)
        angle = VECTNORM(args); axis = UNITVECT(args)
        cos = COS(angle); sin = SIN(angle)
        if axis[1]==axis[2]==0.0:    # rotation about x
            mat[1,1] = cos;    mat[1,2] = -sin;
            mat[2,1] = sin;    mat[2,2] = cos;
        elif axis[0]==axis[2]==0.0:    # rotation about y
            mat[0,0] = cos;    mat[0,2] = sin;
            mat[2,0] = -sin;    mat[2,2] = cos;
        elif axis[0]==axis[1]==0.0:    # rotation about z
            mat[0,0] = cos;    mat[0,1] = -sin;
            mat[1,0] = sin;    mat[1,1] = cos;
        
        else:        # general 3D rotation (Rodrigues' rotation formula)    
            I = scipy.identity(3) ; u = axis
            Ux = scipy.array([
                [0,        -u[2],      u[1]],
                [u[2],        0,     -u[0]],
                [-u[1],     u[0],         0]])
            UU = scipy.array([
                [u[0]*u[0],    u[0]*u[1],    u[0]*u[2]],
                [u[1]*u[0],    u[1]*u[1],    u[1]*u[2]],
                [u[2]*u[0],    u[2]*u[1],    u[2]*u[2]]])
            mat[:3,:3] = cos*I + sin*Ux + (1.0-cos)*UU
        
    
    return mat.view(Mat)

def larEmbed(k):
    def larEmbed0(model):
        V,CV = model
        if k>0:
            V = [v+[0.]*k for v in V] 
        elif k<0:
            V = [v[:-k] for v in V] 
        return V,CV
    return larEmbed0

def larApply(affineMatrix):
    def larApply0(model):
        if isinstance(model,Model):
            # V = scipy.dot([v.tolist()+[1.0] for v in model.verts], affineMatrix.T).tolist()
            V = scipy.dot(array([v+[1.0] for v in model.verts]), affineMatrix.T).tolist()
            V = [v[:-1] for v in V]
            CV = copy.copy(model.cells)
            return Model((V,CV))
        elif isinstance(model,tuple) or isinstance(model,list):
            V,CV = model
            V = scipy.dot([v+[1.0] for v in V], affineMatrix.T).tolist()
            return [v[:-1] for v in V],CV
    return larApply0

""" Flatten a list using Python generators """
def flatten(lst):
    for x in lst:
        if (isinstance(x,tuple) or isinstance(x,list)) and len(x)==2:
            yield x
        elif (isinstance(x,tuple) or isinstance(x,list)):
            for x in flatten(x):
                yield x
        elif isinstance(x, Struct):
            for x in flatten(x.body):
                yield x
        else:
            yield x
 
#  lst = [[1], 2, [[3,4], 5], [[[]]], [[[6]]], 7, 8, []]
#  print list(flatten(lst)) 
#  [1, 2, 3, 4, 5, 6, 7, 8]

#  import itertools
#  chain = itertools.chain.from_iterable([[1,2],[3],[5,89],[],[6]])
#  print(list(chain))
#  [1, 2, 3, 5, 89, 6]    ###  TODO: Bug coi dati sopra?

def checkStruct(lst):
    """ Return the common dimension of structure elements.

        TODO: aggiungere test sulla dimensione minima delle celle (legata a quella di immersione)
    """
    vertsDims = [computeDim(obj) for obj in flatten(lst)]
    vertsDims = [dim for dim in vertsDims if dim!=None and dim!=0]
    if EQ(vertsDims) and len(vertsDims)!=0: 
        return vertsDims[0]
    else: 
        print "\n vertsDims =", vertsDims
        print "*** LAR ERROR: Struct dimension mismatch."

def computeDim(obj):
    """ Check for dimension of a structure element (Verts or V). 
    """
    if (isinstance(obj,Model)):
        return obj.n
    elif (isinstance(obj,tuple) or isinstance(obj,list)) and len(obj)==2:
        V = obj[0]
        if (isinstance(V,list) and isinstance(V[0],list) and 
                (isinstance(V[0][0],float) or isinstance(V[0][0],int))): 
            dim = len(V[0])
            return dim
    elif (isinstance(obj,Mat)):
        dim = obj.shape[0]-1
        return dim
    else: return 0

""" Traversal of a scene multigraph """
def traversal(CTM, stack, obj, scene=[]):
    for i in range(len(obj)):
        if isinstance(obj[i],Model): 
            scene += [larApply(CTM)(obj[i])]
        elif (isinstance(obj[i],tuple) or isinstance(obj[i],list)) and len(obj[i])==2:
            scene += [larApply(CTM)(obj[i])]
        elif isinstance(obj[i],Mat): 
            CTM = scipy.dot(CTM, obj[i])
        elif isinstance(obj[i],Struct):
            stack.append(CTM) 
            traversal(CTM, stack, obj[i], scene)
            CTM = stack.pop()
    return scene

def evalStruct(struct):
    dim = checkStruct(struct.body)
    CTM, stack = scipy.identity(dim+1), []
    scene = traversal(CTM, stack, struct, []) 
    return scene

""" class definitions for LAR """
import scipy
class Mat(scipy.ndarray): pass
class Verts(scipy.ndarray): pass

class Model:
    """ A pair (geometry, topology) of the LAR package """
    def __init__(self,(verts,cells)):
        self.n = len(verts[0])
        # self.verts = scipy.array(verts).view(Verts)
        self.verts = verts
        self.cells = cells
    def __getitem__(self,i):
        return list((self.verts,self.cells))[i]

from myfont import *
class Struct:
    """ The assembly type of the LAR package """
    def __init__(self,data,name=None):
        self.body = data
        if name != None: 
            self.name = str(name)
        else:
            self.name = str(id(self))
        self.box = box(self) 
    def __name__(self):
        return self.name
    def __iter__(self):
        return iter(self.body)
    def __len__(self):
        return len(list(self.body))
    def __getitem__(self,i):
        return list(self.body)[i]
    def __print__(self): 
        return "<Struct name: %s>" % self.__name__()
    def __repr__(self):
        return "<Struct name: %s>" % self.__name__()
        #return "'Struct(%s,%s)'" % (str(self.body),str(str(self.__name__())))
    def draw(self,color=WHITE,scaling=1):
      vmin,vmax = self.box
      delta = VECTDIFF([vmax,vmin])
      point = CCOMB(self.box)
      scalingFactor = scaling*delta[0]/20.
      text = TEXTWITHATTRIBUTES (TEXTALIGNMENT='centre', TEXTANGLE=0,
               TEXTWIDTH=0.1*scalingFactor, 
               TEXTHEIGHT=0.2*scalingFactor,
               TEXTSPACING=0.025*scalingFactor)
      return T([1,2,3])(point)(COLOR(color)(text(self.name)))

""" Structure to pair (Vertices,Cells) conversion """

def struct2lar(structure):
    listOfModels = evalStruct(structure)
    vertDict = dict()
    index,defaultValue,CW,W = -1,-1,[],[]
        
    for model in listOfModels:
        if isinstance(model,Model):
            V,FV = model.verts,model.cells
        elif (isinstance(model,tuple) or isinstance(model,list)) and len(model)==2:
            V,FV = model
        for k,incell in enumerate(FV):
            outcell = []
            for v in incell:
                key = vcode(V[v])
                if vertDict.get(key,defaultValue) == defaultValue:
                    index += 1
                    vertDict[key] = index
                    outcell += [index]
                    W += [eval(key)]
                else: 
                    outcell += [vertDict[key]]
            CW += [outcell]
            
    return W,CW

def larEmbed(k):
    def larEmbed0(model):
        V,CV = model
        if k>0:
            V = [v+[0.]*k for v in V] 
        elif k<0:
            V = [v[:-k] for v in V] 
        return V,CV
    return larEmbed0

""" Computation of the containment box of a Lar Struct or Model """
import copy
def box(model):
    if isinstance(model,Mat): return []
    elif isinstance(model,Struct):
        dummyModel = copy.deepcopy(model)
        dummyModel.body = [term if (not isinstance(term,Struct)) else [term.box,[[0,1]]]  for term in model.body]
        listOfModels = evalStruct( dummyModel )
        dim = len(listOfModels[0][0][0])
        theMin,theMax = box(listOfModels[0]) 
        for theModel in listOfModels[1:]:
            modelMin, modelMax = box(theModel)
            theMin = [val if val<theMin[k] else theMin[k] for k,val in enumerate(modelMin)]
            theMax = [val if val>theMax[k] else theMax[k] for k,val in enumerate(modelMax)]
        return [theMin,theMax]
    elif isinstance(model,Model):
        V = model.verts
    elif (isinstance(model,tuple) or isinstance(model,list)) and len(model)==2:
        V = model[0]
    coords = TRANS(V)
    theMin = [min(coord) for coord in coords]
    theMax = [max(coord) for coord in coords]
    return [theMin,theMax]

