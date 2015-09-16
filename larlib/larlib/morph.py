""" LAR implementation of morphological operators on multidimensional images."""
""" Initial import of modules """
from larlib import *
import scipy.misc, numpy, pickle
from numpy.random import randint


def randomImage(shape, structure, noiseFraction=0.0):
   """ Generation of random image of given shape and structure. 
      Return scipy.ndarray(shape)
   """
   print "noiseFraction =",noiseFraction
   ranges = [shape[k]/structure[k] for k in range(len(shape))]
   random_array = randint(0, 255, size=structure)
   image_array = numpy.zeros(shape)
   for index in scipy.array(CART(AA(range)(shape))):
      block = index/ranges
      if random_array[tuple(block)] < 127:
         image_array[tuple(index)] = 0 
      else: 
         image_array[tuple(index)] = 255
   
   noiseQuantity = PROD(list(shape))*noiseFraction
   k = 0
   while k < noiseQuantity:
      index = tuple(AA(randint)(list(shape)))
      if image_array[index] == 0: image_array[index] = 255
      else: image_array[index] = 0
      k += 1
   if len(shape)==3:
      for k in range(shape[0]):
         scipy.misc.imsave('tmp/outfile'+str(k).zfill(3)+'.png', image_array[k])
   else:
      scipy.misc.imsave('tmp/outfile'+str(k).zfill(3)+'.png', image_array)
   
   return image_array

def dump(object,filename):
   with open(filename, 'wb') as f:
       pickle.dump(object, f)
       
def load(filename):
   with open(filename, 'rb') as f:
       object = pickle.load(f)
   return object

def imageChainComplex (shape):
   tokens = str(shape)[1:-1].split(',')
   tokens = [token.strip() for token in tokens]
   filename = "tmp/larimage-" + "-".join(tokens) + ".pickle"
   
   if os.path.isfile(filename):
      shape, skeletons, operators = loadImageLAR(filename)
   else:
      skeletons = gridSkeletons(list(shape))
      operators = boundaryOps(skeletons)
      imageLAR = (shape, skeletons, operators)
      dump(imageLAR,filename)
      print "filename =",filename
   return shape, skeletons, operators
   
def loadImageLAR(filename):
   object = load(filename)
   shape, skeletons, operators = object
   return shape, skeletons, operators

def mapTupleToInt(shape):
   d = len(shape)
   weights = [PROD(shape[(k+1):]) for k in range(d-1)]+[1]
   
   def mapTupleToInt0(tuple):
      return INNERPROD([tuple,weights])
   return mapTupleToInt0

def setMaskWindow(window,image_array):
   minPoint, maxPoint = window
   imageShape = list(image_array.shape)
   indexRanges = zip(minPoint,maxPoint)
   tuples = CART([range(min,max) for min,max in indexRanges])
   
   imageCochain = image_array.reshape(PROD(imageShape))
   mapping = mapTupleToInt(imageShape)
   windowChain = [mapping(tuple) for tuple in tuples]
   segmentChain = [cell for cell in windowChain if imageCochain[cell]==255]
   
   for cell in segmentChain: imageCochain[cell] = 127
   image_array = imageCochain.reshape(imageShape)
   #for k in range(shape[0]):
   #  scipy.misc.imsave('tmp/outfile'+str(k).zfill(3)+'.png', image_array[k])
   
   return segmentChain

def larImage(shape):
   """ Compute vertices and skeletons of an image of given shape """
   imageVerts = larImageVerts(shape)
   skeletons = gridSkeletons(list(shape))
   return imageVerts, skeletons

def boundaryOps(skeletons):
   """ CSR matrices of boundary operators from list of skeletons """
   return [boundary(skeletons[k+1],faces) 
      for k,faces in enumerate(skeletons[:-1])]


def chainTransform(skeletons,d,chain,csrMatrix,n):
   # n: number of vertices shared in the incidence relation
   cellNumber = len(skeletons[d])
   csrChain = scipy.sparse.csr_matrix((cellNumber,1))
   for h in chain: csrChain[h,0] = 1
   cooOutChain = matrixProduct(csrMatrix, csrChain).tocoo()
   outChain = [cooOutChain.row[h]
      for h,val in enumerate(cooOutChain.data) if int(val) >= n]
   return outChain 

def visImageChain (shape,chain, imageVerts, skeletons):
   # imageVerts, skeletons = larImage(shape)
   chainLAR = [cell for k,cell in enumerate(skeletons[-1]) if k in chain]
   return imageVerts,chainLAR

def imageChainBoundary(shape, operators):
   imageVerts, skeletons = larImage(shape)
   # operators = boundaryOps(skeletons)
   cellNumber = PROD(list(shape))
   
   def imageChainBoundary0(k):
      csrBoundaryMat = operators[-1]
      facets = skeletons[k-1]
      
      def imageChainBoundary1(chain):
         csrChain = scipy.sparse.csr_matrix((cellNumber,1))
         for h in chain: csrChain[h,0] = 1
         csrBoundaryChain = matrixProduct(csrBoundaryMat, csrChain)
         for h,value in enumerate(csrBoundaryChain.data):
            if MOD([value,2]) == 0: csrBoundaryChain.data[h] = 0
         cooBoundaryChain = csrBoundaryChain.tocoo()
         boundaryChain = [cooBoundaryChain.row[h] 
            for h,val in enumerate(cooBoundaryChain.data) if val == 1]
         
         boundaryChainModel = imageVerts, [facets[h] for h in boundaryChain]     
         return boundaryChainModel,boundaryChain
      
      return imageChainBoundary1
   return imageChainBoundary0

def larDown(skeletons,d):
   """ Down operator, to multiply a d-chain and return the incident (d-1)-chain """
   csrMd = csrCreate(skeletons[d])
   csrMinus = csrCreate(skeletons[d-1])
   csrDown = matrixProduct(csrMinus,csrTranspose(csrMd))
   return csrDown
   
def larUp(skeletons,d):
   """ Up operator, to multiply a d-chain and return the incident (d+1)-chain """
   csrMd = csrCreate(skeletons[d])
   csrPlus = csrCreate(skeletons[d+1])
   csrUp = matrixProduct(csrPlus,csrTranspose(csrMd))
   return csrUp   

def UUD(skeletons,d):
   """ Compute the morphological operator UUP.
      Return a CSR matrix to be applied to the coordinate representation of a chain
   """ 
   D = larDown(skeletons,d)
   U1 = larUp(skeletons,d-1)
   U2 = larUp(skeletons,d)
   UUDout = matrixProduct(U2,matrixProduct(U1,D))
   return D,U1,U2,UUDout

def testAlgebraicMorphologyStepByStep (solid, b_rep, chain, imageVerts, skeletons):
   d = len(skeletons)-1
   D,U1,U2,csrMorphOp = UUD(skeletons,d-1)
   
   VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(solid)))
   VIEW(COLOR(MAGENTA)(STRUCT(MKPOLS(b_rep))))
   
   outputChain = chainTransform(skeletons,d-1,chain,D,d-1)
   chainLAR = [cell for k,cell in enumerate(skeletons[0]) if k in outputChain]
   model0 = (imageVerts,chainLAR)
   VIEW(COLOR(RED)(STRUCT(MKPOLS(model0))))
   
   M2 = csrCreate(skeletons[d])
   
   outputChain = chainTransform(skeletons,0,outputChain,M2,d)
   chainLAR = [cell for k,cell in enumerate(skeletons[d]) if k in outputChain]
   model2 = (imageVerts,chainLAR)
   VIEW(COLOR(YELLOW)(STRUCT(MKPOLS(model2))))
   
   return imageVerts,chainLAR
 
def testAlgebraicMorphology (solid, b_rep, chain, imageVerts, skeletons):
   d = len(skeletons)-1
   D,U1,U2,csrMorphOp = UUD(skeletons,d-1)
   outputChain = chainTransform(skeletons,1,chain,csrMorphOp,d)
   chainLAR = [cell for k,cell in enumerate(skeletons[d]) if k in outputChain]
   model2 = (imageVerts,chainLAR)
   return imageVerts,chainLAR
 
