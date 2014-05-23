\documentclass[11pt,oneside]{article}	%use"amsart"insteadof"article"forAMSLaTeXformat
\usepackage{geometry}		%Seegeometry.pdftolearnthelayoutoptions.Therearelots.
\geometry{letterpaper}		%...ora4paperora5paperor...
%\geometry{landscape}		%Activateforforrotatedpagegeometry
%\usepackage[parfill]{parskip}		%Activatetobeginparagraphswithanemptylineratherthananindent
\usepackage{graphicx}				%Usepdf,png,jpg,orepsßwithpdflatex;useepsinDVImode
								%TeXwillautomaticallyconverteps-->pdfinpdflatex		
\usepackage{amssymb}
\usepackage[colorlinks]{hyperref}

\usepackage{framed}
\usepackage{amsthm}
\newtheorem{remark}{Remark}
\newtheorem{definition}{Definition}

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

\title{Imaging Morphology with LAR
\footnote{This document is part of the \emph{Linear Algebraic Representation with CoChains} (LAR-CC) framework~\cite{cclar-proj:2013:00}. \today}
}
\author{Alberto Paoluzzi}
%\date{}							%Activatetodisplayagivendateornodate

\begin{document}
\maketitle
\nonstopmode

\begin{abstract}
In this module we aim to implement the four operators of mathematical morphology, i.e.~the \emph{dilation}, \emph{erosion}, \emph{opening} and \emph{closing} operators, by the way of matrix operations representing the linear operators---\emph{boundary} and \emph{coboundary}---over LAR. 
According to the multidimensional character of LAR, our implementation is dimension-independent.
In few words, it works as follows: (a)  the input is (the coordinate representation of) a $d$-chain $\gamma$; (b) compute its boundary $\partial_d(\gamma)$; (c) extract the maximal $(d-2)$-chain $\epsilon \subset \partial_d(\gamma)$; (d) consider the $(d-1)$-chain returned from its coboundary $\delta_{d-2}(\epsilon)$; (e) compute the $d$-chain $\eta := \delta_{d-1}(\delta_{d-2}(\epsilon)) \subset C_d$ \emph{without} performing the  $\mbox{mod\ 2}$ final transformation on the resulting coordinate vector, that would provide a zero result, according to the standard algebraic constraint $\delta\circ\delta=0$. It is easy to show that $\eta \equiv (\oplus \gamma) - (\ominus \gamma)$ provides the \emph{morphological gradient} operator. The four standard morphological operators are therefore  consequently computable.
\end{abstract}

\tableofcontents

\section{Test image generation}

Various methods for the input or the generation of a test image  are developed in the subsections of this section. The aim is to prepare a set of controlled test beds, used to check both the implementation and the working properties of our topological implementation of morphological operators. 


\subsection{Random binary multidimensional image}

A multidimensional binary image is generated here by using a random approach, both for the bulk structure and the small artefacts of the image.  


%-------------------------------------------------------------------------------
@d Generation of random image
@{def randomImage(shape, structure, noiseFraction=0.0):
	""" Generation of random image of given shape and structure. 
		Return scipy.ndarray(shape)
	"""
	print "noiseFraction =",noiseFraction
	@< Generation of bulk array structure @>
	@< Generation of random artifacts @>
	return image_array
@}
%-------------------------------------------------------------------------------


\paragraph{Generation of the gross image}
First we generate a 2D grid of squares by Cartesian product, and produce the bulk of the random image then used to test our approach to morphological operators via topological ones.


%import scipy
%scipy.misc.imsave('outfile.jpg', image_array)
%
%scipy.ndimage.imread(fname, flatten=False, mode=None)[source]
%
%rand(10, 10)
%
	
%-------------------------------------------------------------------------------
@d Generation of bulk array structure
@{ranges = [shape[k]/structure[k] for k in range(len(shape))]
random_array = randint(0, 255, size=structure)
image_array = numpy.zeros(shape)
for index in scipy.array(CART(AA(range)(shape))):
	block = index/ranges
	if random_array[tuple(block)] < 127:
		image_array[tuple(index)] = 0 
	else: 
		image_array[tuple(index)] = 255
@}
%-------------------------------------------------------------------------------



\paragraph{Generation of random artefacts upon the image}

Then random noise is added to the previously generated image, in order to produce artifacts at the pixel scale. 

%-------------------------------------------------------------------------------
@d Generation of random artifacts
@{noiseQuantity = PROD(list(shape))*noiseFraction
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
@}
%-------------------------------------------------------------------------------


\section{Selection of an image segment}

In this section we implement several methods for image segmentation and segment selection. 

\subsection{Selection of a test chain}

The first and simplest method is the selection of the portion of a binary image contained within a masking window.
Here we select the (white) sub-image contained in a given window, and compute the coordinate representation of the (chain) sub-image.

\paragraph{Mask definition}

A \emph{window} within a $d$-image is defined by $2\times d$ integer numbers (2 multi-indices), corresponding to the window  \texttt{minPoint} (minimum indices) and to the window \texttt{maxPoint} (maximum indices). A list of multi-index tuples, contained in the \texttt{window} variable, is generated by the function \texttt{setMaskWindow} below.

%-------------------------------------------------------------------------------
@d Generation of a masking window
@{def setMaskWindow(window,image_array):
	minPoint, maxPoint = window
	imageShape = list(image_array.shape)
	@< Generation of multi-index window @>
	@< Window-to-chain mapping @>
	@< Change chain color to grey @>
	return segmentChain
@}
%-------------------------------------------------------------------------------

The set of tuples of indices contained in a (multidimensional) window is given below.
 
%-------------------------------------------------------------------------------
@d Generation of multi-index window
@{indexRanges = zip(minPoint,maxPoint)
tuples = CART([range(min,max) for min,max in indexRanges])
@}
%-------------------------------------------------------------------------------



\subsection{Mapping of integer tuples to integers}

In order to produce the coordinate representation of a chain in a multidimensional image (or $d$-image) we need: (a) to choose a basis of image elements, i.e.~of $d$-cells, and in particular to fix an ordering of them; (b) to map the multidimensional index, selecting a single $d$-cell of the image, to a single integer mapping the cell to its linear position within the chosen basis ordering. 


\paragraph{Grid of hyper-cubes of unit size}
Let $S_i=(0,1,...,n_i-1)$ be ordered integer sets with $n_i$ elements, and 
\[
S= S_0 \times S_1 \times \cdots \times S_{d-1}
\] 
the set of indices of elements of a $d$-image.

\begin{definition}[$d$-image shape]
The \emph{shape} of a $d$-image with $n_0\times n_1 \times\cdots\times n_{d-1}$ elements (here called \emph{voxels}) is the ordered set $(n_0, n_1, \ldots, n_{d-1})$.
\end{definition}


\paragraph{$d$-dimensional row-major order}

Given a $d$-image with shape $S=(n_0,n_1,...,n_{d-1})$ and number of elements $n=\prod n_i$, 
the mapping  
\[
S_0 \times S_1 \times \cdots \times S_{d-1} \to \{ 0, 1, \ldots, n-1\}
\]
 is a {linear combination} with integer {weights}  $(w_0,w_1,...,w_{d-2},1)$, such that:
\[
(i_0,i_1,...,i_{d-1}) \mapsto i_0 w_0 +i_1 w_1 +\cdots +i_{d-1} w_{d-1},
\]
where 
\[
w_k = n_{k+1}  n_{k+2} \cdots  n_{d-1}, \qquad 0\leq k\leq d-2.
\]

\paragraph{Implementation}
A functional implementation of the \emph{Tuples to integers mapping} is given by the second-order  \texttt{mapTupleToInt} function, that  accepts in a first application the \texttt{shape} of the image (to compute the tuple space of indices of $d$-cells), and then takes a single tuple in the second application. Of course, the function  returns the cell address in the linear address space associated to the given \texttt{shape}.

%-------------------------------------------------------------------------------
@d Tuples to integers mapping
@{def mapTupleToInt(shape):
	d = len(shape)
	weights = [PROD(shape[(k+1):]) for k in range(d-1)]+[1]
	
	def mapTupleToInt0(tuple):
		return INNERPROD([tuple,weights])
	return mapTupleToInt0
@}
%-------------------------------------------------------------------------------


\paragraph{From tuples multi-indices to chain coordinates}

The set of address \texttt{tuples} of $d-cells$ ($d$-dimensional image elements) within the \emph{mask} is here mapped to the corresponding set of (single) integers associated to the low-level image elements (pixels or voxels, depending on the image dimension and shape), denoted \texttt{windowChain}. Such total chain of the mask \texttt{window} is then filtered to contain the only coordinates of \emph{white} image elements within the window, and returned as the set of integer cell indices \texttt{segmentChain}.


%-------------------------------------------------------------------------------
@d Window-to-chain mapping
@{imageCochain = image_array.reshape(PROD(imageShape))
mapping = mapTupleToInt(imageShape)
windowChain = [mapping(tuple) for tuple in tuples]
segmentChain = [cell for cell in windowChain if imageCochain[cell]==255]
@}
%-------------------------------------------------------------------------------

\subsection{Show segment chain on binary image}

Now we need to show visually the selected \texttt{segmentChain}, by change the color of its cells from white (255) to middle grey (127). Just remember that \texttt{imageCochain} is the linear representation of the image, with number of cells equal to \texttt{PROD(imageShape)}. Then the modified image is restored within \texttt{image\_array}, and is finally exported to a \texttt{.png} image file.

%-------------------------------------------------------------------------------
@d Change chain color to grey
@{for cell in segmentChain: imageCochain[cell] = 127
image_array = imageCochain.reshape(imageShape)
#for k in range(shape[0]):
#	scipy.misc.imsave('tmp/outfile'+str(k).zfill(3)+'.png', image_array[k])
@}
%-------------------------------------------------------------------------------


\section{Construction of (co)boundary operators}

A $d$-image is a \emph{cellular $d$-complex} where cells are $k$-cuboids ($0\leq k\leq d$), i.e.~Cartesian products of a number $k$ of 1D intervals, embedded in $d$-dimensional Euclidean space. 


\subsection{Reading and writing LAR of image from disk}
Since the construction of the \emph{chain complex} supported by a cellular complex with $O(n^d)$ $d$-cells---and $n=O(10^3)$---may be really time-consuming, it is unquestionably useful to store on disk, once and for all, the topological model of a multidimensional image, i.e.~the LAR of its chain complex,  and just restore it (or the needed part of it --- TODO) when necessary, since \emph{it only depends on the \texttt{shape} of a $d$-image}, i.e.~from the array arrangement of its $d$-cells (``hypervoxels'').

%-------------------------------------------------------------------------------
@d Save and restore data object from file
@{def dump(object,filename):
	with open(filename, 'wb') as f:
	    pickle.dump(object, f)
	    
def load(filename):
	with open(filename, 'rb') as f:
	    object = pickle.load(f)
	return object
@}
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
@d Save or restore the chain complex of multidimensional image
@{def imageChainComplex (shape):
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
@}
%-------------------------------------------------------------------------------
	



\subsection{LAR chain complex construction}

In our first multidimensional implementation of morphological operators through algebraic topology of the image seen as a cellular complex, we compute the whole sequence of characteristic matrices $M_k$ ($0\leq k \leq d$) in \texttt{BRC} form, and the whole sequence of matrices $[\partial_k]$ ($0\leq k \leq d$) in \texttt{CSR} form.

\begin{figure}[htbp] %  figure placement: here, top, bottom, or page
   \centering
   \includegraphics[width=0.7\linewidth]{images/larcomplex} 
   \caption{The LAR definition of  a chain complex: a sequence of characteristic matrices \emph{and} a sequence of boundary operators.}
   \label{fig:example}
\end{figure}

\paragraph{Array of characteristic matrices}

A direct construction of cuboidal complexes is offered, within the \texttt{larcc} package, by the \texttt{largrid.larCuboids} function. 

%-------------------------------------------------------------------------------
@d Characteristic matrices of multidimensional image
@{def larImage(shape):
	""" Compute vertices and skeletons of an image of given shape """
	imageVerts = larImageVerts(shape)
	skeletons = gridSkeletons(list(shape))
	return imageVerts, skeletons
@}
%-------------------------------------------------------------------------------

\paragraph{Example}
Consider a (very!) small 3D image of \texttt{shape=(2,2,2)}. The data structures returned by the \texttt{larImage} function are shown below, where \texttt{imageVerts} gives the integer coordinates of vertices of the 3D (image) complex, and \texttt{skeletons} is the list of characteristic matrices $M_k$ ($0\leq k\leq d$) in \texttt{BRC} form.

%-------------------------------------------------------------------------------
@D Example of characteristic matrices (and vertices) of multidimensional image
@{imageVerts, skeletons = larImage((2,2,2))

print imageVerts,
>>> [[0,0,0],[0,0,1],[0,0,2],[0,1,0],[0,1,1],[0,1,2],[0,2,0],[0,2,1],[0,2,2],
[1,0,0],[1,0,1],[1,0,2],[1,1,0],[1,1,1],[1,1,2],[1,2,0],[1,2,1],[1,2,2],[2,0,
0],[2,0,1],[2,0,2],[2,1,0],[2,1,1],[2,1,2],[2,2,0],[2,2,1],[2,2,2]]

print skeletons[1:],
>>> [
[[0,1],[1,2],[3,4],[4,5],[6,7],[7,8],[9,10],[10,11],[12,13],[13,14],[15,
16],[16,17],[18,19],[19,20],[21,22],[22,23],[24,25],[25,26],[0,3],[1,4],[2,
5],[3,6],[4,7],[5,8],[9,12],[10,13],[11,14],[12,15],[13,16],[14,17],[18,21],
[19,22],[20,23],[21,24],[22,25],[23,26],[0,9],[1,10],[2,11],[3,12],[4,13],[5,
14],[6,15],[7,16],[8,17],[9,18],[10,19],[11,20],[12,21],[13,22],[14,23],[15,
24],[16,25],[17,26]],
[[0,1,3,4],[1,2,4,5],[3,4,6,7],[4,5,7,8],[9,10,12,13],[10,11,13,14],[12,13,15,
16],[13,14,16,17],[18,19,21,22],[19,20,22,23],[21,22,24,25],[22,23,25,26],[0,
1,9,10],[1,2,10,11],[3,4,12,13],[4,5,13,14],[6,7,15,16],[7,8,16,17],[9,10,18,
19],[10,11,19,20],[12,13,21,22],[13,14,22,23],[15,16,24,25],[16,17,25,26],[0,
3,9,12],[1,4,10,13],[2,5,11,14],[3,6,12,15],[4,7,13,16],[5,8,14,17],[9,12,18,
21],[10,13,19,22],[11,14,20,23],[12,15,21,24],[13,16,22,25],[14,17,23,26]],
[[0,1,3,4,9,10,12,13],[1,2,4,5,10,11,13,14],[3,4,6,7,12,13,15,16],[4,5,7,8,13,
14,16,17],[9,10,12,13,18,19,21,22],[10,11,13,14,19,20,22,23],[12,13,15,16,21,
22,24,25],[13,14,16,17,22,23,25,26]]
]
@}
%-------------------------------------------------------------------------------


\paragraph{Array of matrices of boundary operators}

The function \texttt{boundaryOps} takes the array of \texttt{BRC} reprs of characteristic matrices, and returns the array of \texttt{CSR} matrix reprs of boundary operators $\partial_k$ ($1\leq k\leq d$).

%-------------------------------------------------------------------------------
@d CSR matrices of boundary operators
@{def boundaryOps(skeletons):
	""" CSR matrices of boundary operators from list of skeletons """
	return [boundary(skeletons[k+1],faces) 
		for k,faces in enumerate(skeletons[:-1])]
@}
%-------------------------------------------------------------------------------


\paragraph{Boundary chain of a $k$-chain of a $d$-image}

%-------------------------------------------------------------------------------
@d Boundary of image chain computation
@{def imageChainBoundary(shape, operators):
	imageVerts, skeletons = larImage(shape)
	# operators = boundaryOps(skeletons)
	cellNumber = PROD(list(shape))
	
	def imageChainBoundary0(k):
		csrBoundaryMat = operators[-1]
		facets = skeletons[k-1]
		
		def imageChainBoundary1(chain):
			@< Boundary*chain product and interpretation @>
			boundaryChainModel = imageVerts, [facets[h] for h in boundaryChain]		
			return boundaryChainModel,boundaryChain
		
		return imageChainBoundary1
	return imageChainBoundary0
@}
%-------------------------------------------------------------------------------

\paragraph{Low-level SpMSpV matrix-chain product}

A low-level implementation of the product of boundary matrix times the coordinate representation of the image chain under consideration is given below. It was enveloped and protected in this macro because it is strongly dependent on the CSR structures and tools provided by the sparse matrix module of the \texttt{scipy} library.	

%-------------------------------------------------------------------------------
@d Boundary*chain product and interpretation
@{csrChain = scipy.sparse.csr_matrix((cellNumber,1))
for h in chain: csrChain[h,0] = 1
csrBoundaryChain = matrixProduct(csrBoundaryMat, csrChain)
for h,value in enumerate(csrBoundaryChain.data):
	if MOD([value,2]) == 0: csrBoundaryChain.data[h] = 0
cooBoundaryChain = csrBoundaryChain.tocoo()
boundaryChain = [cooBoundaryChain.row[h] 
	for h,val in enumerate(cooBoundaryChain.data) if val == 1]
@}
%-------------------------------------------------------------------------------


\subsection{Visualisation of an image chain and its boundary}


\paragraph{$d$-Chain visualisation}

The \texttt{visImageChain} function given by the macro \emph{Visualisation of an image chain} below. 

%-------------------------------------------------------------------------------
@d Pyplasm visualisation of an image chain
@{def visImageChain (shape,chain, imageVerts, skeletons):
	# imageVerts, skeletons = larImage(shape)
	chainLAR = [cell for k,cell in enumerate(skeletons[-1]) if k in chain]
	return imageVerts,chainLAR
@}
%-------------------------------------------------------------------------------


\begin{figure}[htbp] %  figure placement: here, top, bottom, or page
   \centering
   \includegraphics[height=0.3\linewidth,width=0.16\linewidth]{images/morph-chain2D1} 
   \includegraphics[height=0.3\linewidth,width=0.16\linewidth]{images/morph-chain1D1} 
   \includegraphics[height=0.3\linewidth,width=0.16\linewidth]{images/morph-chain2D2} 
   \includegraphics[height=0.3\linewidth,width=0.16\linewidth]{images/morph-chain1D2} 
   \includegraphics[height=0.3\linewidth,width=0.16\linewidth]{images/morph-chain2D3} 
   \includegraphics[height=0.3\linewidth,width=0.16\linewidth]{images/morph-chain1D3} 
   \caption{example caption}
   \label{fig:morph}
\end{figure}


\section{The bulk of morphological imaging}

The main part of the implementation of our algebraic construction of morphological operators on $d$-images is given in this section. First we introduce two generalized incidence operators, respectively named \emph{$\mathcal{U}_d$}, which stands for \emph{Upper}, and \emph{$\mathcal{L}_d$}, which stands for \emph{Lower}. Then we implement the standard morphological operators of Dilation (\emph{$\mathcal{D}_d$}), Erosion (\emph{$\mathcal{E}_d$}), Opening (\emph{$\mathcal{O}_d$}), and Closing (\emph{$\mathcal{C}_d$}).

\subsection{``Down'' and ``Up'' operators on Chains}

The Up and Down incidence operators are defined as mapping $d$-cells to $(d+1)$- and $(d+1)$-cells that, respectively, share vertices with them:
\[
\mathcal{U}_d : C_d \to C_{d+1}, \qquad \mathcal{D}_d : C_d \to C_{d-1}.
\]

Remembering that the characteristic matrix $M_d$ corresponds to the mapping from vertices to $d$-cells, we have
\[
[\mathcal{U}_d] = M_{d+1}\,M_d^t, \qquad [\mathcal{D}_d] = M_{d-1}\,M_d^t.
\]

\paragraph{Down and Up operators} Just consider that the characteristic matrices $M_{d-1},M_{d},M_{d+1}$,  stored in position $d-1$, $d$, and $d+1$ of the \texttt{skeletons} array, respectively provide the subsets of vertices of proper dimensions incising on each cell. They are used to compute the topological incidence operators, according to the paper~\cite{Dicarlo:2014:TNL:2543138.2543294}.

%-------------------------------------------------------------------------------
@d Down and Up operators
@{def larDown(skeletons,d):
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
@}
%-------------------------------------------------------------------------------

\paragraph{The $\mathcal{UUD}$ operator}
A generic $\mathcal{UUD}: C_{d-1} \to C_{d}$ operator is given to be associated with the boundary $\partial_2$, in order to return the $\mathcal{ERO}\cup\mathcal{DIL}$ (TODO:correct) of a $d$-chain.

%-------------------------------------------------------------------------------
@d UUD operator: maps Down-UP-UP its input chains
@{def UUD(skeletons,d):
	""" Compute the morphological operator UUP.
		Return a CSR matrix to be applied to the coordinate representation of a chain
	""" 
	D = larDown(skeletons,d)
	U1 = larUp(skeletons,d-1)
	U2 = larUp(skeletons,d)
	UUDout = matrixProduct(U2,matrixProduct(U1,D))
	return D,U1,U2,UUDout
@}
%-------------------------------------------------------------------------------

A pair of similar function \texttt{testAlgebraicMorphology} and \texttt{testAlgebraicMorphologyStepByStep} are given here, both accepting the same (redundant) inputs, in order to test our algebraic approach to mathematical morphology. In particular: 

%\item
%\texttt{solid} is the LAR model (i.e.~the pair (\texttt{V},\texttt{CV})) of the original image;
%\item
%\texttt{b_rep} is the LAR model of its boundary (\texttt{V},\texttt{BV});
%\item
%\texttt{chain} is the $d$-chain of the image, given as the list of indices of $(d-1)$-cells within the boundary's $(d-1)$-chain;
%\item
%\texttt{imageVerts} 
%\item
%\texttt{skeletons} 

%-------------------------------------------------------------------------------
@d Testing the algebraic morphology
@{def testAlgebraicMorphology (solid, b_rep, chain, imageVerts, skeletons):
	d = len(skeletons)-1
	D,U1,U2,csrMorphOp = UUD(skeletons,d-1)
	outputChain = chainTransform(skeletons,1,chain,csrMorphOp,d)
	chainLAR = [cell for k,cell in enumerate(skeletons[d]) if k in outputChain]
	model2 = (imageVerts,chainLAR)
	return imageVerts,chainLAR
@}
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
@d Testing the algebraic morphology method step-by-step
@{def testAlgebraicMorphologyStepByStep (solid, b_rep, chain, imageVerts, skeletons):
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
@}
%-------------------------------------------------------------------------------

\paragraph{Boundary chain mapping}
The function \texttt{chainTransform} takes as input a \texttt{chain} of dimension \texttt{d} (actually the boundary of a chain of dimension $d+1$), the \texttt{skeletons} array of \texttt{BRC} representations of matrices $M_k$ ($0\leq k\leq d+1$), and the number $n$ that $(d+1)$-cells resulting by multiplication with the operator's matrix \texttt{csrMatrix} must ahare with they incident $d$-cells in order to be included in the output of the considered algebraic morphology operator.

%-------------------------------------------------------------------------------
@d csrMatrix*chain product and interpretation
@{
def chainTransform(skeletons,d,chain,csrMatrix,n):
	# n: number of vertices shared in the incidence relation
	cellNumber = len(skeletons[d])
	csrChain = scipy.sparse.csr_matrix((cellNumber,1))
	for h in chain: csrChain[h,0] = 1
	cooOutChain = matrixProduct(csrMatrix, csrChain).tocoo()
	outChain = [cooOutChain.row[h]
		for h,val in enumerate(cooOutChain.data) if int(val) >= n]
	return outChain 
@}
%-------------------------------------------------------------------------------


\paragraph{Example}

%-------------------------------------------------------------------------------
@d Algebraic dilation minus erosion
@{D1 = larDown(skeletons,1)
U0 = larUp(skeletons,0)
U1 = larUp(skeletons,1)
UUD1 = UUD(skeletons,1)
@}
%-------------------------------------------------------------------------------



\subsection{Maximal ($d-2$)-chain extraction}


%-------------------------------------------------------------------------------
@d Extract the maximal ($d-2$)-chain from a ($d-1$)-chain
@{

@}
%-------------------------------------------------------------------------------

\subsection{($d-1$)-Star of a ($d-2$)-chain computation}

%-------------------------------------------------------------------------------
@d Compute the ($d-1$)-star of a ($d-2$)-chain
@{

@}
%-------------------------------------------------------------------------------


\subsection{$d$-Star of a $(d-1)$-chain computation}

%-------------------------------------------------------------------------------
@d Compute the d-star of a ($d-1$)-chain
@{

@}
%-------------------------------------------------------------------------------


\subsection{Dilation and erosion computation}



\subsection{Opening and closing computation}




\section{Exporting the \texttt{morph} module}

\paragraph{Exporting the morph module}
%-------------------------------------------------------------------------------
@o lib/py/morph.py
@{""" LAR implementation of morphological operators on multidimensional images."""
@< Initial import of modules @>
@< Generation of random image @>
@< Save and restore data object from file @>
@< Save or restore the chain complex of multidimensional image @>
@< Tuples to integers mapping @>
@< Generation of a masking window @>
@< Characteristic matrices of multidimensional image @>
@< CSR matrices of boundary operators @>
@< csrMatrix*chain product and interpretation @>
@< Pyplasm visualisation of an image chain @>
@< Boundary of image chain computation @>
@< Down and Up operators @>
@< UUD operator: maps Down-UP-UP its input chains @>
@< Testing the algebraic morphology method step-by-step @> 
@< Testing the algebraic morphology @> 
@}
%-------------------------------------------------------------------------------

The set of importing commends needed by test files in this module is given in the macro below.

%-------------------------------------------------------------------------------
@d Initial import of modules
@{import sys,os
import scipy.misc, numpy, pickle
from numpy.random import randint
from pyplasm import *

""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')

@< Import the module @(largrid@) @>
@}
%-------------------------------------------------------------------------------


\section{Morphological operations examples}

\subsection{2D image masking and boundary computation}

\paragraph{Test example}

The \texttt{larcc.morph} API is used here to generate a random black and white image, with an \emph{image segment} selected and extracted by masking, then colored in middle grey, and exported to an image file.  
The \texttt{shape} and \texttt{structure} variables must contain two tuples, of equal \texttt{len}, of powers of 2, used toy define the size of random blocks in the generated image. The \texttt{window} variable is used to define the portion of the image where the \emph{segmentChain} of white pixels (or voxels, for 3-images) is computed.

%------------------------------------------------------------------
@o test/py/morph/test01.py
@{@< Initial import of modules @>
@< Import the module @(morph@) @> 
shape = 64,64
structure = 8,8
assert len(shape) == len(structure)
imageVerts = larImageVerts(shape)
_, skeletons, operators = imageChainComplex (shape)
image_array = randomImage(shape, structure, 0.05)
minPoint, maxPoint = (0,0), (64,64)
window = minPoint, maxPoint
segmentChain = setMaskWindow(window,image_array)
	
solid = visImageChain (shape,segmentChain, imageVerts, skeletons)
b_rep,boundaryChain = imageChainBoundary(shape, operators)(2)(segmentChain)

stepwiseTest = True
if stepwiseTest:
	model = testAlgebraicMorphologyStepByStep(solid, b_rep, 
				boundaryChain, imageVerts, skeletons)
else:
	model = testAlgebraicMorphology(solid, b_rep, 
				boundaryChain, imageVerts, skeletons)
VIEW(STRUCT(MKPOLS(model)))

V = model[0]
M = AA(tuple)(model[1])
S = AA(tuple)(solid[1])
B = AA(tuple)(b_rep[1])
D = list(set(S).union(M))
E = list(set(S).difference(M))
M,S,B,D,E = (AA(AA(list)))([M,S,B,D,E])
M2 = STRUCT(MKPOLS((V,M)))
S2 = STRUCT(MKPOLS((V,S)))

S2 = COLOR(CYAN)(STRUCT(MKPOLS((V,S))))
B1 = COLOR(MAGENTA)(STRUCT(MKPOLS((V,B))))
D2 = COLOR(YELLOW)(STRUCT(MKPOLS((V,D))))
E2 = COLOR(WHITE)(STRUCT(MKPOLS((V,E))))
VIEW(STRUCT([D2,S2,E2,B1]))

VIEW(STRUCT([D2,S2,B1]))
VIEW(STRUCT([S2,E2,B1]))

@}
%------------------------------------------------------------------

%A0 = matrix([
%[0,2,0],
%[4,0,0],
%[0,0,2],
%])
%B0 = matrix([
%[0,1,0],
%[1,0,0],
%[0,0,1],
%])
%A1 = matrix([
%[0,10,0],
%[0,11,0],
%[5,0,12],
%])
%B1 = matrix([
%[0,1,0],
%[0,1,0],
%[1,0,1],
%])


%===============================================================================
\appendix
\section{Utilities}

\subsection{Importing a generic module}
First we define a parametric macro to allow the importing of \texttt{larcc} modules from the project repository \texttt{lib/py/}. When the user needs to import some project's module, she may call this macro as done in Section~\ref{sec:lar2psm}.
%------------------------------------------------------------------
@d Import the module
@{import @1
from @1 import *
@}
%------------------------------------------------------------------

\paragraph{Importing a module} A function used to import a generic \texttt{lacccc} module within the current environment is also useful.
%------------------------------------------------------------------
@d Function to import a generic module
@{def importModule(moduleName):
	@< Import the module @(moduleName@) @>
@| importModule @}
%------------------------------------------------------------------


\bibliographystyle{amsalpha}
\bibliography{morph}

\end{document}
