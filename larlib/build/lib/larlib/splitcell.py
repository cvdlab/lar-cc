from pyplasm import *

# COVECTOR
def COVECTOR(points):
  pointdim=len(points[0])
  plane=Planef.bestFittingPlane(pointdim,[item for sublist in points for item in sublist])
  return [plane.get(I) for I in range(0,pointdim+1)]

# SPLITCELL
def SPLITCELL(plane,points,tolerance=1e-6,ntry=4):
  navigator=GraphNavigator()
  pointdim=len(points[0])
  Vmat,Hmat=[Matf(),Matf()]
  graph=Graph.mkpol(Vmat,Hmat,pointdim,len(points),[item for sublist in points for item in sublist],tolerance)
  graph.transform(Vmat,Hmat)
  #assert(graph.getNCells(2)==1)
  [graph,cell]=[graph,graph.each(pointdim).getNode()]
  [below,equal,above]=graph.split(navigator,cell,Planef(plane),tolerance,ntry)
  assert(below>=0 and equal>=0 and above>=0) #otherwise failed
  
  def findPoints(cell) :
    N=graph.findCells(0,cell,navigator) 
    points=[graph.getVecf(navigator.getCell(0,I)) for I in range(0,N)]
    return [[p.get(I) for I in range(1,pointdim+1)] for p in points]  
    
  return [findPoints(below),findPoints(equal),findPoints(above)]

if __name__ == "__main__":
    plane = COVECTOR([[0.0,0.5],  [1.0,0.5],  [2.0,0.5],  [3.0,0.5]]) # example of plane Y=0.5
    plane = COVECTOR([[0.0,0.0],  [1.0,0.0],  [2.0,0.0],  [3.0,0.0]]) # example of plane Y=0.0
    plane = COVECTOR([[0.0,1.0],  [1.0,1.0],  [2.0,1.0],  [3.0,1.0]]) # example of plane Y=0.1
    plane = COVECTOR([[0.0,2.0],  [1.0,2.0],  [2.0,2.0],  [3.0,2.0]]) # example of plane Y=0.2
    cell  = [[0.0,0.0],  [1.0,0.0],  [1.0,1.0],  [0.0,1.0]] # example of cuboid in range [0,1]^2
    [below,equal,above]=SPLITCELL(plane,cell)
    print "plane",plane
    print "below",below
    print "equal",equal
    print "above",above

