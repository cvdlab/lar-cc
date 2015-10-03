"""Module with functions needed to interface LAR with pyplasm via Struct"""
from larlib import *

""" TODO: use package Decimal (http://docs.python.org/2/library/decimal.html) """
global PRECISION
PRECISION = 4.

def verySmall(number): return abs(number) < 10**-(PRECISION)

def prepKey (args): return "["+", ".join(args)+"]"

def fixedPrec(value):
    out = round(value*10**(PRECISION))/10**(PRECISION)
    if out == -0.0: out = 0.0
    return str(out)
    
def vcode (vect): 
    """
    To generate a string representation of a number array.
    Used to generate the vertex keys in PointSet dictionary, and other similar operations.
    """
    return prepKey(AA(fixedPrec)(vect))

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
        if len(model)==2: V,CV = model
        elif len(model)==3: V,CV,FV = model
        if k>0:
            V = [v+[0.]*k for v in V] 
        elif k<0:
            V = [v[:-k] for v in V] 
        if len(model)==2: return V,CV
        elif len(model)==3: return V,CV,FV
    return larEmbed0

def larEmbed(k):
    def larEmbed0(model):
        if k>0:
            model[0] = [v+[0.]*k for v in model[0]] 
        elif k<0:
            model[0] = [v[:-k] for v in model[0]] 
        return model
    return larEmbed0

""" Apply an affine transformation to a LAR model  """
from scipy import array

def larApply(affineMatrix):
    def larApply0(model):
        if isinstance(model,Model):
            V = scipy.dot(array([v+[1.0] for v in model.verts]), affineMatrix.T).tolist()
            V = [v[:-1] for v in V]
            CV = copy.copy(model.cells)
            return Model((V,CV))
        elif isinstance(model,tuple) or isinstance(model,list):
            if len(model)==2: V,CV = model
            elif len(model)==3: V,CV,FV = model
            V = scipy.dot([list(v)+[1.0] for v in V], affineMatrix.T).tolist()
            if len(model)==2: return [v[:-1] for v in V],CV
            elif len(model)==3: return [v[:-1] for v in V],CV,FV
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
    obj = lst[0]
    if (isinstance(obj,tuple) or isinstance(obj,list)):
        dim = len(obj[0][0])
    elif isinstance(obj,Model): 
        dim = obj.n    
    elif isinstance(obj,Mat): 
        dim = obj.shape[0]-1    
    elif isinstance(obj,Struct): 
        dim = len(obj.box[0])    
    return dim

""" Remove duplicate faces  """
from collections import defaultdict
def removeDups (CW):
    CW = list(set(AA(tuple)(CW)))
    CWs = list(set(AA(tuple)  (AA(sorted)(CW))  ))
    no_duplicates = defaultdict(list)
    for f in CWs: no_duplicates[f] = []
    for f in CW:
        no_duplicates[tuple(sorted(f))] += [f]
    CW = [f[0] for f in no_duplicates.values()]
    return CW

""" Remove the unused vertices """
def larRemoveVertices(V,FV):
    vertDict = dict()
    index,defaultValue,FW,W = -1,-1,[],[]
        
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
        FW += [outcell]
    return W,FW

""" Traversal of a scene multigraph """
def traversal(CTM, stack, obj, scene=[]):
    for i in range(len(obj)):
        if isinstance(obj[i],Model): 
            scene += [larApply(CTM)(obj[i])]
        elif (isinstance(obj[i],tuple) or isinstance(obj[i],list)) and (
                len(obj[i])==2 or len(obj[i])==3):
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

""" Struct iterable class """
class Struct:
    """ The assembly type of the LAR package """
    def __init__(self,data=None,name=None,category=None):
        if data==None or data==[]:
            self.body = []
        else:
            self.body = [item for item in data if item != None]
            self.box = box(self) 
            self.dim = len(self.box[0])
        if name != None: 
            self.name = str(name)
        else:
            self.name = str(id(self))
        if category != None: 
            self.category = str(category)
        else:
            self.category = "feature"
    def __name__(self):
        return self.name
    def __category__(self):
        return self.category
    def __iter__(self):
        return iter(self.body)
    def __len__(self):
        return len(list(self.body))
    def __getitem__(self,i):
        return list(self.body)[i]
    def __setitem__(self,i,value):
        self.body[i] = value
    def __print__(self): 
        return "<Struct name: %s>" % self.__name__()
    def __repr__(self):
        return "<Struct name: %s>" % self.__name__()
        #return "'Struct(%s,%s)'" % (str(self.body),str(str(self.__name__())))
    def set_name(self,name):
        self.name = str(name)
    def clone(self,i=0):
        from copy import deepcopy
        newObj = deepcopy(self)
        if i != 0: newObj.name = self.name + "_" + str(i)
        return newObj
    def set_category(self,category):
        self.category = str(category)
    def boundary(self):
        data = struct2lar(self)
        if len(data) == 3:
            V,FV,EV = data
            #import pdb; pdb.set_trace()
            return V,FV,EV
        else:
            return "<Struct name: %s> boundary non computable" % self.__name__()
    def draw(self,color=WHITE,scaling=1,metric=ID):
        vmin,vmax = self.box
        delta = VECTDIFF([vmax,vmin])
        point = CCOMB(self.box)
        scalingFactor = scaling*delta[0]/20.
        text = TEXTWITHATTRIBUTES (TEXTALIGNMENT='centre', TEXTANGLE=0,
                    TEXTWIDTH=0.1*scalingFactor, 
                    TEXTHEIGHT=0.2*scalingFactor,
                    TEXTSPACING=0.025*scalingFactor)
        point = metric([point])[0]
        return T([1,2,3])(point)(COLOR(color)(text(self.name)))

""" Structure to pair (Vertices,Cells) conversion """

def struct2lar(structure,metric=ID):
    listOfModels = evalStruct(structure)
    vertDict = dict()
    index,defaultValue,CW,W,FW = -1,-1,[],[],[]
        
    for model in listOfModels:
        if isinstance(model,Model):
            V,FV = model.verts,model.cells
        elif (isinstance(model,tuple) or isinstance(model,list)):
            if len(model)==2: V,FV = model
            elif len(model)==3: V,FV,EV = model
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
        if len(model)==3:
            for k,incell in enumerate(EV):
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
                FW += [outcell]
            
    if ((isinstance(model,tuple) or isinstance(model,list)) and len(model)==2) or (
        (isinstance(model,Model) and model.n==2)): 
        if len(CW[0])==2: 
            CW = list(set(AA(tuple)(AA(sorted)(CW))))
        else: CW = removeDups(CW)
        return metric(W),CW
    if ((isinstance(model,tuple) or isinstance(model,list)) and len(model)==3) or (
        (isinstance(model,Model) and model.n==3)): 
        FW = list(set(AA(tuple)(AA(sorted)(FW))))
        CW = removeDups(CW)
        return metric(W),CW,FW

def larEmbed(k):
    def larEmbed0(model):
        if len(model)==2: V,CV = model
        elif len(model)==3: V,CV,FV = model
        if k>0:
            V = [v+[0.]*k for v in V] 
        elif k<0:
            V = [v[:-k] for v in V] 
        if len(model)==2: return V,CV
        elif len(model)==3: return V,CV,FV
    return larEmbed0

def larEmbed(k):
    def larEmbed0(model):
        if k>0:
            model[0] = [v+[0.]*k for v in model[0]] 
        elif k<0:
            model[0] = [v[:-k] for v in model[0]] 
        return model
    return larEmbed0

""" embed a struct object """
""" Structure embedding algorithm """
def embedTraversal(cloned, obj,n,suffix):
    for i in range(len(obj)):
        if isinstance(obj[i],Model): 
            cloned.body += [obj[i]]
        elif (isinstance(obj[i],tuple) or isinstance(obj[i],list)) and (
                len(obj[i])==2):
            V,EV = obj[i]
            V = [v+n*[0.0] for v in V]
            cloned.body  += [(V,EV)]
        elif (isinstance(obj[i],tuple) or isinstance(obj[i],list)) and (
                len(obj[i])==3):
            V,FV,EV = obj[i]
            V = [v+n*[0.0] for v in V]
            cloned.body  += [(V,FV,EV)]
        elif isinstance(obj[i],Mat): 
            mat = obj[i]
            d,d = mat.shape

            newMat = scipy.identity(d+n*1)
            for h in range(d-1): 
                for k in range(d-1): 
                    newMat[h,k] = mat[h,k]
                newMat[h,d-1+n*1] = mat[h,d-1]
            cloned.body  +=  [newMat.view(Mat)]

        elif isinstance(obj[i],Struct):
            newObj = Struct()
            newObj.box = hstack((obj[i].box, [n*[0],n*[0]]))
            newObj.name = obj[i].name+suffix
            newObj.category = obj[i].category
            cloned.body  += [embedTraversal(newObj, obj[i], n, suffix)]
    return cloned


def embedStruct(n):
    def embedStruct0(struct,suffix="New"):
        if n==0: 
            return struct, len(struct.box[0])
        cloned = Struct()
        cloned.box = hstack((struct.box, [n*[0],n*[0]])).tolist()
        cloned.name = str(id(cloned))  #struct.name+suffix
        cloned.category = struct.category
        cloned.dim = struct.dim + n
        cloned = embedTraversal(cloned,struct,n,suffix) 
        return cloned
    return embedStruct0

""" Computation of the containment box of a Lar Struct or Model """
import copy
def box(model):
    if isinstance(model,Mat): return []
    elif isinstance(model,Struct):
        dummyModel = copy.deepcopy(model)
        dummyModel.body = [term if (not isinstance(term,Struct)) else [term.box,[[0,1]]]  for term in model.body]
        listOfModels = evalStruct( dummyModel )
        #dim = checkStruct(listOfModels)
        theMin,theMax = box(listOfModels[0]) 
        for theModel in listOfModels[1:]:
            modelMin, modelMax = box(theModel)
            theMin = [val if val<theMin[k] else theMin[k] for k,val in enumerate(modelMin)]
            theMax = [val if val>theMax[k] else theMax[k] for k,val in enumerate(modelMax)]
        return [theMin,theMax]
    elif isinstance(model,Model):
        V = model.verts
    elif (isinstance(model,tuple) or isinstance(model,list)) and (len(model)==2 or len(model)==3):
        V = model[0]
    coords = TRANS(V)
    theMin = [min(coord) for coord in coords]
    theMax = [max(coord) for coord in coords]
    return [theMin,theMax]

