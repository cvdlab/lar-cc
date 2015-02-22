\documentclass[11pt,oneside]{article}    %use"amsart"insteadof"article"forAMSLaTeXformat
\usepackage{geometry}        %Seegeometry.pdftolearnthelayoutoptions.Therearelots.
\geometry{letterpaper}        %...ora4paperora5paperor...
%\geometry{landscape}        %Activateforforrotatedpagegeometry
%\usepackage[parfill]{parskip}        %Activatetobeginparagraphswithanemptylineratherthananindent
\usepackage{graphicx}                %Usepdf,png,jpg,orepsßwithpdflatex;useepsinDVImode
                                %TeXwillautomaticallyconverteps-->pdfinpdflatex        
\usepackage{amssymb}
\usepackage[colorlinks]{hyperref}

%----macros begin---------------------------------------------------------------
\usepackage{color}
\usepackage{amsthm}

\def\conv{\mbox{\textrm{conv}\,}}
\def\aff{\mbox{\textrm{aff}\,}}
\def\E{\mathbb{E}}
\def\R{\mathbb{R}}
\def\Z{\mathbb{Z}}
\def\tex{\TeX}
\def\latex{\LaTeX}
\def\v#1{{\bf #1}}
\def\p#1{{\bf #1}}
\def\T#1{{\bf #1}}

\def\vet#1{{\left(\begin{array}{cccccccccccccccccccc}#1\end{array}\right)}}
\def\mat#1{{\left(\begin{array}{cccccccccccccccccccc}#1\end{array}\right)}}

\def\lin{\mbox{\rm lin}\,}
\def\aff{\mbox{\rm aff}\,}
\def\pos{\mbox{\rm pos}\,}
\def\cone{\mbox{\rm cone}\,}
\def\conv{\mbox{\rm conv}\,}
\newcommand{\homog}[0]{\mbox{\rm homog}\,}
\newcommand{\relint}[0]{\mbox{\rm relint}\,}

%----macros end-----------------------------------------------------------------

\title{Boolean chains
\footnote{This document is part of the \emph{Linear Algebraic Representation with CoChains} (LAR-CC) framework~\cite{cclar-proj:2013:00}. \today}
}
\author{Alberto Paoluzzi}
%\date{}                            %Activatetodisplayagivendateornodate

\begin{document}
\maketitle
\nonstopmode

\begin{abstract}
A novel algorithm for computation of Boolean operations between cellular complexes is given in this module.
It is based on bucketing of possibly interacting geometry using a box-extension of kd-trees, normally used for  point proximity queries. 
Such kd-tree representation of containment boxes of cells, allow us to compute a number of independent buckets of data,  with some elements possibly replicated,  to be used for brute force intersection, followed by elimination of duplicated data.
We use a fast algorithm~\cite{Bajaj:1996:SCC:237218.237246} for splitting a cell with an hyperplane, and append the splitted sub-cells to their buckets, in order to apply all the splits induced by the boundary cells contained in the same bucket. 
A final tagging of cells as either belonging or not to each operand follows, allowing for fast extraction of Boolean results between any pair of chains (subsets of cells).
This Boolean algorithm can be considered of a \emph{Map-Reduce} kind, and hence suitable of a distributed implementation over big datasets. The actual engineered implementation will follow the present prototype, using some distributed NoSQL database, like MongoDB or Riak. 
\end{abstract}

\tableofcontents

%===============================================================================
\section{Introduction}
%===============================================================================


%===============================================================================
\section{Preview of the algorithm}
%===============================================================================

The whole Boolean algorithm is composed by four stages in sequence, denoted in the following as \emph{Unification}, \emph{Bucketing}, \emph{Intersection}, and \emph{Reconstruction}. The algorithm described here is both multidimensional and variadic. Multidimensional means that the arguments are solid in Euclidean space of dimension $d$, with $d$ small integer.
The \emph{arity}  of a function or operation is the number of arguments or operands the function or operation accepts. 
In computer science, a function accepting a variable number of arguments is called \emph{variadic}.

\subsection{Unification}
%===============================================================================

In this first step the boundaries of the $n$ Boolean arguments are computed and merged together as a set of chains defined in the discrete set $\texttt{V}$ made by the union of their vertices, and possibly by a set of random points. 

Random points are only added when the $L_1$ norm (Manhattan or taxicab norm) of $(d-1)$-dimensional boundary cells is greater than a given fraction of the bounding box $BB(<args>)$ of the vertices of arguments.

The Delaunay triangulation $\texttt{CV}$ is hence computed, as well as its $(d-1)$ skeleton \texttt{FV}, providing the  \texttt{LAR} model \texttt{V,CV} of an initial \emph{partition} of the space. 

It is important to notice here that the  cell decompositions of the $n$ arguments, and in particula their boundaries, have vertices belonging to the same discrete space. The arguments of the boolean expression provide a \emph{covering} of the same space.

Actually, only the (\emph{oriented}) boundaries \texttt{V,FV$_i$} $(1\leq i\leq n)$ of the varius arguments are retained, togheter with the \texttt{V,CV} triangulation, for the successive steps of the algorithm.

\subsection{Bucketing}
%===============================================================================

The bounding boxes of \texttt{CV} and \texttt{FV$_i$} are computed, and their \emph{box-kd-tree} is worked-out, so providing a group of buckets of close cells, that can be elaborated independently, and possibly in parallel, to compute the intersection of the Delaunay simplices with the boundary cells.  

\subsection{Intersection}
%===============================================================================

\subsection{Reconstruction}
%===============================================================================



%===============================================================================
\section{Implementation}
%===============================================================================

\subsection{Box-kd-tree}
%===============================================================================


\paragraph{Split the boxes between the below,above subsets}
%-------------------------------------------------------------------------------
@D Split the boxes between the below,above subsets
@{""" Split the boxes between the below,above subsets """
def splitOnThreshold(boxes,subset,coord):
    theBoxes = [boxes[k] for k in subset]
    threshold = centroid(theBoxes,coord)
    ncoords = len(boxes[0])/2
    a = coord%ncoords
    b = a+ncoords
    below,above = [],[]
    for k in subset:
        if boxes[k][a] <= threshold: below += [k]
    for k in subset:
        if boxes[k][b] >= threshold: above += [k]
    return below,above
@}
%-------------------------------------------------------------------------------

\paragraph{Test if bucket OK or append to splitting stack}
%-------------------------------------------------------------------------------
@D Test if bucket OK or append to splitting stack
@{""" Test if bucket OK or append to splitting stack """
def splitting(bucket,below,above, finalBuckets,splittingStack):
    if (len(below)<4 and len(above)<4) or len(set(bucket).difference(below))<7 \
        or len(set(bucket).difference(above))<7: 
        finalBuckets.append(below)
        finalBuckets.append(above)
    else: 
        splittingStack.append(below)
        splittingStack.append(above)
@}
%-------------------------------------------------------------------------------


\paragraph{Remove subsets from bucket list}
%-------------------------------------------------------------------------------
@D Remove subsets from bucket list @{
""" Remove subsets from bucket list """
def removeSubsets(buckets):
    n = len(buckets)
    A = zeros((n,n))
    for i,bucket in enumerate(buckets):
        for j,bucket1 in enumerate(buckets):
            if set(bucket).issubset(set(bucket1)):
                A[i,j] = 1
    B = AA(sum)(A.tolist())
    out = [bucket for i,bucket in enumerate(buckets) if B[i]==1]
    return out

def geomPartitionate(boxes,buckets):
    geomInters = [set() for h in range(len(boxes))]
    for bucket in buckets:
        for k in bucket:
            geomInters[k] = geomInters[k].union(bucket)
    for h,inters in enumerate(geomInters):
        geomInters[h] = geomInters[h].difference([h])
    return AA(list)(geomInters)
@}
%-------------------------------------------------------------------------------
    


\paragraph{Iterate the splitting until \texttt{splittingStack} is empty}
%-------------------------------------------------------------------------------
@D Iterate the splitting until splittingStack is empty
@{""" Iterate the splitting until \texttt{splittingStack} is empty """
def boxTest(boxes,h,k):
    B1,B2,B3,B4,B5,B6,_ = boxes[k]
    b1,b2,b3,b4,b5,b6,_ = boxes[h]
    return not (b4<B1 or B4<b1 or b5<B2 or B5<b2 or b6<B3 or B6<b3)

def boxBuckets(boxes):
    bucket = range(len(boxes))
    splittingStack = [bucket]
    finalBuckets = []
    while splittingStack != []:
        bucket = splittingStack.pop()
        below,above = splitOnThreshold(boxes,bucket,1)
        below1,above1 = splitOnThreshold(boxes,above,2)
        below2,above2 = splitOnThreshold(boxes,below,2) 
               
        below11,above11 = splitOnThreshold(boxes,above1,3)
        below21,above21 = splitOnThreshold(boxes,below1,3)        
        below12,above12 = splitOnThreshold(boxes,above2,3)
        below22,above22 = splitOnThreshold(boxes,below2,3)  
              
        splitting(above1,below11,above11, finalBuckets,splittingStack)
        splitting(below1,below21,above21, finalBuckets,splittingStack)
        splitting(above2,below12,above12, finalBuckets,splittingStack)
        splitting(below2,below22,above22, finalBuckets,splittingStack)
        
        finalBuckets = list(set(AA(tuple)(finalBuckets)))
    parts = geomPartitionate(boxes,finalBuckets)
    parts = [[h for h in part if boxTest(boxes,h,k)] for k,part in enumerate(parts)]
    return AA(sorted)(parts)
@}
%-------------------------------------------------------------------------------

\paragraph{aaaaaa}
%-------------------------------------------------------------------------------
@D aaaaaa
@{""" aaaaa """

@}
%-------------------------------------------------------------------------------


\subsection{Merging the boundaries}
%===============================================================================

\subsection{Elementary splitting}
%===============================================================================

\subsection{Boolean chains}
%===============================================================================

%===============================================================================
\section{Esporting the Library}
%===============================================================================

%-------------------------------------------------------------------------------
@O lib/py/bool2.py
@{""" Module for Boolean computations between geometric objects """
from pyplasm import *
""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')
from inters import *
DEBUG = True

@< Coding utilities @>
@< Split the boxes between the below,above subsets @>
@< Test if bucket OK or append to splitting stack @>
@< Remove subsets from bucket list @>
@< Iterate the splitting until splittingStack is empty @>
@}
%-------------------------------------------------------------------------------

%===============================================================================
\section{Test examples}
%===============================================================================

\subsection{Random triangles}
%===============================================================================


\paragraph{Generation of random triangles and their boxes}
%-------------------------------------------------------------------------------
@O test/py/bool2/test01.py
@{""" Generation of random triangles and their boxes """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *
glass = MATERIAL([1,0,0,0.1,  0,1,0,0.1,  0,0,1,0.1, 0,0,0,0.1, 100])

randomTriaArray = randomTriangles(10,0.99)
VIEW(STRUCT(AA(MKPOL)([[verts, [[1,2,3]], None] for verts in randomTriaArray])))

boxes = containmentBoxes(randomTriaArray)
hexas = AA(box2exa)(boxes)
cyan = COLOR(CYAN)(STRUCT(AA(MKPOL)([[verts, [[1,2,3]], None] for verts in randomTriaArray])))
yellow = STRUCT(AA(glass)(AA(MKPOL)([hex for hex,qualifier in hexas])))
VIEW(STRUCT([cyan,yellow]))
@}
%-------------------------------------------------------------------------------


\paragraph{Generation of random quadrilaterals and their boxes}
%-------------------------------------------------------------------------------
@O test/py/bool2/test02.py
@{""" Generation of random quadrilaterals and their boxes """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *
glass = MATERIAL([1,0,0,0.1,  0,1,0,0.1,  0,0,1,0.1, 0,0,0,0.1, 100])

randomQuadArray = randomQuads(10,1)
VIEW(STRUCT(AA(MKPOL)([[verts, [[1,2,3,4]], None] for verts in randomQuadArray])))

boxes = containmentBoxes(randomQuadArray)
hexas = AA(box2exa)(boxes)
cyan = COLOR(CYAN)(STRUCT(AA(MKPOL)([[verts, [[1,2,3,4]], None] for verts in randomQuadArray])))
yellow = STRUCT(AA(glass)(AA(MKPOL)([hex for hex,qualifier in hexas])))
VIEW(STRUCT([cyan,yellow]))
@}
%-------------------------------------------------------------------------------




\subsection{Testing the box-kd-trees}
%===============================================================================


\paragraph{Visualizing with different colors the buckets of box-kd-tree}
%-------------------------------------------------------------------------------
@O test/py/bool2/test04.py @{
""" Visualizing with different colors the buckets of box-kd-tree """
from pyplasm import *
""" import modules from larcc/lib """
import sys
sys.path.insert(0, 'lib/py/')
from bool2 import *
glass = MATERIAL([1,0,0,0.1,  0,1,0,0.1,  0,0,1,0.1, 0,0,0,0.1, 100])

randomQuadArray = randomQuads(30,0.8)
VIEW(STRUCT(AA(MKPOL)([[verts, [[1,2,3,4]], None] for verts in randomQuadArray])))

boxes = containmentBoxes(randomQuadArray)
hexas = AA(box2exa)(boxes)
glass = MATERIAL([1,0,0,0.1,  0,1,0,0.1,  0,0,1,0.1, 0,0,0,0.1, 100])
yellow = STRUCT(AA(glass)(AA(MKPOL)([hex for hex,data in hexas])))
VIEW(STRUCT([#cyan,
    yellow]))

parts = boxBuckets(boxes)
for k,part in enumerate(parts):
    bunch = [glass(STRUCT( [MKPOL(hexas[h][0]) for h in part]))]
    bunch += [COLOR(RED)(MKPOL(hexas[k][0]))]
    VIEW(STRUCT(bunch))
@}
%-------------------------------------------------------------------------------

bunch=[]
for k,part in enumerate(parts):
    bunch += [[h for h in part if boxTest(boxes,h,k)]]

\appendix
%===============================================================================
\section{Code utilities}
%===============================================================================

\paragraph{Coding utilities}

Some utility fuctions used by the module are collected in this appendix. Their macro names can be seen in the below script.

%-------------------------------------------------------------------------------
@D Coding utilities
@{""" Coding utilities """
@< Generation of a random 3D point @>
@< Generation of random 3D triangles @>
@< Generation of random 3D quadrilaterals @>
@< Generation of a single random triangle @>
@< Containment boxes @>
@< Transformation of a 3D box into an hexahedron @>
@< Computation of the 1D centroid of a list of 3D boxes @>
@}
%-------------------------------------------------------------------------------


\paragraph{Generation of random triangles}
The function \texttt{randomTriangles} returns the array \texttt{randomTriaArray} with a given number of triangles generated within the unit 3D interval. The \texttt{scaling} parameter is used to scale every such triangle, generated by three randow points, that could be possibly located to far from each other, even at the distance of the diagonal of the unit cube.

The arrays \texttt{xs}, \texttt{ys} and \texttt{zs}, that contain the $x,y,z$ coordinates of triangle points, are used to compute the minimal translation \texttt{v} needed to transport the entire set of data within the positive octant of the 3D space. 

%-------------------------------------------------------------------------------
@D Generation of random 3D triangles
@{""" Generation of random triangles """
def randomTriangles(numberOfTriangles=400,scaling=0.3):
    randomTriaArray = [rtriangle(scaling) for k in range(numberOfTriangles)]
    [xs,ys,zs] = TRANS(CAT(randomTriaArray))
    xmin, ymin, zmin = min(xs), min(ys), min(zs)
    v = array([-xmin,-ymin, -zmin])
    randomTriaArray = [[list(v1+v), list(v2+v), list(v3+v)] for v1,v2,v3 in randomTriaArray]
    return randomTriaArray
@}
%-------------------------------------------------------------------------------

\paragraph{Generation of random 3D quadrilaterals}

%-------------------------------------------------------------------------------
@D Generation of random 3D quadrilaterals
@{""" Generation of random 3D quadrilaterals """
def randomQuads(numberOfQuads=400,scaling=0.3):
    randomTriaArray = [rtriangle(scaling) for k in range(numberOfQuads)]
    [xs,ys,zs] = TRANS(CAT(randomTriaArray))
    xmin, ymin, zmin = min(xs), min(ys), min(zs)
    v = array([-xmin,-ymin, -zmin])
    randomQuadArray = [AA(list)([ v1+v, v2+v, v3+v, v+v2-v1+v3 ]) for v1,v2,v3 in randomTriaArray]
    return randomQuadArray
@}
%-------------------------------------------------------------------------------


\paragraph{Generation of a random 3D point}
A single random point, codified in floating point format, and with a fixed (quite small) number of digits, is returned by the \texttt{rpoint()} function, with no input parameters.
%-------------------------------------------------------------------------------
@D Generation of a random 3D point
@{""" Generation of a random 3D point """
def rpoint():
    return eval( vcode([ random.random(), random.random(), random.random() ]) )
@}
%-------------------------------------------------------------------------------
    
\paragraph{Generation of a single random triangle}
A single random triangle, scaled about its centroid by the \texttt{scaling} parameter, is returned by the \texttt{rtriangle()} function, as a tuple ot two random points in the unit square.
%-------------------------------------------------------------------------------
@D Generation of a single random triangle
@{""" Generation of a single random triangle """
def rtriangle(scaling):
    v1,v2,v3 = array(rpoint()), array(rpoint()), array(rpoint())
    c = (v1+v2+v3)/3
    pos = rpoint()
    v1 = (v1-c)*scaling + pos
    v2 = (v2-c)*scaling + pos
    v3 = (v3-c)*scaling + pos
    return tuple(eval(vcode(v1))), tuple(eval(vcode(v2))), tuple(eval(vcode(v3)))
@}
%-------------------------------------------------------------------------------
    

\paragraph{Containment boxes}

Given as input a list \texttt{randomTriaArray} of pairs of 2D points, the function \texttt{containmentBoxes} returns, in the same order, the list of \emph{containment boxes} of the input lines. A \emph{containment box} of a geometric object of dimension $d$ is defined as the minimal $d$-cuboid, equioriented with the reference frame, that contains the object. For a 2D line it is given by the tuple $(x1,y1,x2,y2)$, where $(x1,y1)$ is the point of minimal coordinates, and $(x2,y2)$ is the point of maximal  coordinates.

%-------------------------------------------------------------------------------
@D Containment boxes
@{""" Containment boxes """
def containmentBoxes(randomPointArray,qualifier=0):
    if len(randomPointArray[0])==2:
        boxes = [eval(vcode([min(x1,x2), min(y1,y2), min(z1,z2), 
                             max(x1,x2), max(y1,y2), max(z1,z2)]))+[qualifier]
                for ((x1,y1,z1),(x2,y2,z2)) in randomPointArray]
    elif len(randomPointArray[0])==3:
        boxes = [eval(vcode([min(x1,x2,x3), min(y1,y2,y3), min(z1,z2,z3), 
                             max(x1,x2,x3), max(y1,y2,y3), max(z1,z2,z3)]))+[qualifier]
                for ((x1,y1,z1),(x2,y2,z2),(x3,y3,z3)) in randomPointArray]
    elif len(randomPointArray[0])==4:
        boxes = [eval(vcode([min(x1,x2,x3,x4), min(y1,y2,y3,y4), min(z1,z2,z3,z4), 
                             max(x1,x2,x3,x4), max(y1,y2,y3,y4), max(z1,z2,z3,z4)]))+[qualifier]
                for ((x1,y1,z1),(x2,y2,z2),(x3,y3,z3),(x4,y4,z4)) in randomPointArray]
    return boxes
@}
%-------------------------------------------------------------------------------

    
\paragraph{Transformation of a 3D box into an hexahedron}
The transformation of a 2D box into a closed rectangular polyline, given as an ordered sequwncw of 2D points, is produced by the function \texttt{box2exa}
%-------------------------------------------------------------------------------
@D Transformation of a 3D box into an hexahedron
@{""" Transformation of a 3D box into an hexahedron """    
def box2exa(box):
    x1,y1,z1,x2,y2,z2,type = box
    verts = [[x1,y1,z1], [x1,y1,z2], [x1,y2,z1], [x1,y2,z2], [x2,y1,z1], [x2,y1,z2], [x2,y2,z1], [x2,y2,z2]]
    cell = [range(1,len(verts)+1)]
    return [verts,cell,None],type

def lar2boxes(model,qualifier=0):
    V,CV = model
    boxes = []
    for k,cell in enumerate(CV):
        verts = [V[v] for v in cell]
        x1,y1,z1 = [min(coord) for coord in TRANS(verts)]
        x2,y2,z2 = [max(coord) for coord in TRANS(verts)]
        boxes += [eval(vcode([min(x1,x2),min(y1,y2),min(z1,z2),max(x1,x2),max(y1,y2),max(z1,z2)]))+[(qualifier,k)]]
    return boxes
@}
%-------------------------------------------------------------------------------
    
\paragraph{Computation of the 1D centroid of a list of 3D boxes}
The 1D \texttt{centroid} of a list of 3D boxes is computed by the function given below.
The direction of computation (either $x,y$ or $z$) is chosen depending on the value of the \texttt{coord} parameter. 
%-------------------------------------------------------------------------------
@D Computation of the 1D centroid of a list of 3D boxes
@{""" Computation of the 1D centroid of a list of 3D boxes """    
def centroid(boxes,coord):
    delta,n = 0,len(boxes)
    ncoords = len(boxes[0])/2
    a = coord%ncoords
    b = a+ncoords
    for box in boxes:
        delta += (box[a] + box[b])/2
    return delta/n
@}
%-------------------------------------------------------------------------------

\bibliographystyle{amsalpha}
\bibliography{bool2}

\end{document}
