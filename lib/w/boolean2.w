\documentclass[11pt,oneside]{article}	%use"amsart"insteadof"article"forAMSLaTeXformat
\usepackage{geometry}		%Seegeometry.pdftolearnthelayoutoptions.Therearelots.
\geometry{letterpaper}		%...ora4paperora5paperor...
%\geometry{landscape}		%Activateforforrotatedpagegeometry
%\usepackage[parfill]{parskip}		%Activatetobeginparagraphswithanemptylineratherthananindent
\usepackage{graphicx}				%Usepdf,png,jpg,oreps•withpdflatex;useepsinDVImode
								%TeXwillautomaticallyconverteps-->pdfinpdflatex		
\usepackage{amssymb}
\usepackage{hyperref}

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

\title{Boolean Chains: set operations with topological chains
\footnote{This document is part of the \emph{Linear Algebraic Representation with CoChains} (LAR-CC) framework~\cite{cclar-proj:2013:00}. \today}
}
\author{Alberto Paoluzzi}
%\date{}							%Activatetodisplayagivendateornodate

\begin{document}
\maketitle
%\nonstopmode

%----macros end-----------------------------------------------------------------
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
\begin{abstract}
Boolean operations are a major addition to every geometric package. Union, intersection, difference and complementation of decomposed spaces are discussed and implemented in this module by making use of the Linear Algebraic Representation (LAR) introduced in~\cite{Dicarlo:2014:TNL:2543138.2543294}. First, the two finite decompositions are merged, by merging their vertices (0-cells of support spaces); then a Delaunay complex based on the vertex set union is computed, and the shared $d$-chain is extracted and split, according to the cellular structure of the input $d$-chains. The results of a Boolean operation are finally computed by sum, product or difference of the (binary) coordinate representation of the (split) argument chains, by using the novel chain-basis resulted from the splitting step. Differently from the totality of algorithms known to the authors, neither search nor traversal of some (complicated) data structure is performed by this algorithm. 
\end{abstract}

\tableofcontents
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
\section{Introduction}

In this section we introduce and shortly outline our novel algorithm for Boolean operations with chain of cells from different space decompositions implemented in this LAR-CC software module.

The input objects are denoted in the remainder as $X_1$ and $X_2$, and their finite cell decompositions as $\Lambda^1$ and $\Lambda^1$. Our goal is to compute $X = X_1\, op\, X_2$, where $op \in \{\cup ,\cap , - ,\ominus \}$ or $\complement X$, based on a common decomposition $\Lambda = \Lambda^1\, op\, \Lambda^2$, with $\Lambda$ being a suitably fragmented decomposition of the X space.

Of course, we aim to compute a minimal (in some sense) decomposition, making the best use of the LAR framework, based on CSR representation of sparse binary matrices and standard matrix algebra operations.
However, in this first implementation of the chain approach to Boolean operations, we are satisfied with a solution using simplicial triangulations of input spaces. Future revisions of our algorithm will be based on more general cellular complexes.

\subsection{User interface}

The API will contain the high-level binary functions \texttt{union}, \texttt{intersection}, \texttt{difference}, and \texttt{xor}. Each of them will call the same function \texttt{boolOps} and then suitably operates the two returned bit arrays, i.e.~the coordinate representations of the input spaces in the merged cell decomposition. The input parameters \texttt{lar1} and  \texttt{lar2} stand for two LAR models, each one constituted by a pair \texttt{(V,CV)}, i.e.~by the matrix \texttt{V} of vertex coordinates and by an integer array \texttt{CV} giving the vertex indices of each $d$-cell.

%-------------------------------------------------------------------------------
@D High-level boolean operations
@{def union(lar1,lar2):
	lar = boolOps(lar1,lar2)
def intersection(lar1,lar2):
	lar = boolOps(lar1,lar2)
def difference(lar1,lar2):
	lar = boolOps(lar1,lar2)
def xor(lar1,lar2):
	lar = boolOps(lar1,lar2)
@}
%-------------------------------------------------------------------------------

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
\section{Merging 0-cells}
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

The real work is performed by the function \texttt{boolOps}, that will procede step-by-step (see Section~\ref{sec:splitting}) to the computation of the minimally fragmented common cell complex where to compute the chain resulting from the requested Boolean operation.

%-------------------------------------------------------------------------------
\subsection{Global reordering of vertex coordinates}
%-------------------------------------------------------------------------------
A global reordering of vertex coordinates is executed as the first step of the Boolean algorithm, in order to eliminate the duplicate vertices, by substituting double vertex copies (coming from the two close points) with a single instance. 

Two dictionaries are created, then merged in a single dictionary, and finally split into three subsets of (vertex,index) pairs, with the aim of rebuilding the input representations, by making use of a novel and more useful vertex indexing.

The union set of vertices is finally reordered using the three subsets of vertices belonging (a) only to the first argument, (b) only to the second argument and (c) to both, respectively denoted as $V_1, V_2, V_{12}$. A top-down description of this initial computational step is provided by the set of macros discussed in this section.

%-------------------------------------------------------------------------------
@D Place the vertices of Boolean arguments in a common space
@{""" First step of Boolen Algorithm """
@< Initial indexing of vertex positions @>
@< Merge two dictionaries with keys the point locations @>
@< Filter the common dictionary into three subsets @>
@< Compute an inverted index to reorder the vertices of Boolean arguments @>
@< Return the single reordered pointset and the two $d$-cell arrays @>
@}
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
\subsubsection{Re-indexing of vertices}
%-------------------------------------------------------------------------------

\paragraph{Initial indexing of vertex positions}
The input LAR models are located in a common space by (implicitly) joining \texttt{V1} and \texttt{V2} in a same array, and (explicitly) shifting the vertex indices in \texttt{CV2} by the length of \texttt{V1}.
%-------------------------------------------------------------------------------
@D Initial indexing of vertex positions
@{from collections import defaultdict, OrderedDict

def vertexSieve(model1, model2):
	V1,CV1 = model1; V2,CV2 = model2
	n = len(V1); m = len(V2)
	def shift(CV, n): 
		return [[v+n for v in cell]for cell in CV]
	CV2 = shift(CV2,n)
@}
%-------------------------------------------------------------------------------

\paragraph{Merge two dictionaries with point location as keys}
Since currently \texttt{CV1} and \texttt{CV2} point to a set of vertices larger than their initial sets 
\texttt{V1} and \texttt{V2}, we index the set $\texttt{V1} \cup \texttt{V2}$ using a Python \texttt{defaultdict} dictionary, in order to avoid errors of "missing key". As dictionary keys, we use the string representation of the vertex position vector provided by the \texttt{vcode} function given in the Appendix.
%-------------------------------------------------------------------------------
@D Merge two dictionaries with keys the point locations
@{
	vdict1 = defaultdict(list)
	for k,v in enumerate(V1): vdict1[vcode(v)].append(k) 
	vdict2 = defaultdict(list)
	for k,v in enumerate(V2): vdict2[vcode(v)].append(k+n) 
	
	vertdict = defaultdict(list)
	for point in vdict1.keys(): vertdict[point] += vdict1[point]
	for point in vdict2.keys(): vertdict[point] += vdict2[point]
@}
%-------------------------------------------------------------------------------

\paragraph{Example of string coding of a vertex position}
The position vector of a point of real coordinates is provided by the function \texttt{vcode}.
An example of coding is given below. The \emph{precision} of the string representation can be tuned at will.
{\small
\begin{verbatim}
>>> vcode([-0.011660381062724849, 0.297350056848685860])
'[-0.0116604, 0.2973501]'
\end{verbatim}}



\paragraph{Filter the common dictionary into three subsets}
\texttt{Vertdict}, dictionary of vertices, uses as key stye position vectors of vertices coded as string, and as values the list of integer indices of vertices on the given position. If the point position belongs either to the first or to second argument only, it is stored in \texttt{case1} or \texttt{case2} lists respectively. If the position (\texttt{item.key}) is shared between two vertices, it is stored in \texttt{case12}.
The variables \texttt{n1}, \texttt{n2}, and \texttt{n12} remember the number of vertices respectively stored in each repository.
%-------------------------------------------------------------------------------
@D Filter the common dictionary into three subsets
@{
	case1, case12, case2 = [],[],[]
	for item in vertdict.items():
		key,val = item
		if len(val)==2:  case12 += [item]
		elif val[0] < n: case1 += [item]
		else: case2 += [item]
	n1 = len(case1); n2 = len(case12); n3 = len(case2)
@}
%-------------------------------------------------------------------------------

\paragraph{Compute an inverted index to reorder the vertices of Boolean arguments}
The new indices of vertices are computed according with their position within the storage repositories \texttt{case1}, \texttt{case2}, and \texttt{case12}. Notice that every \texttt{item[1]} stored in \texttt{case1} or \texttt{case2} is a list with only one integer member. Two such values are conversely stored in each \texttt{item[1]} within \texttt{case12}.
%-------------------------------------------------------------------------------
@D Compute an inverted index to reorder the vertices of Boolean arguments 
@{
	invertedindex = list(0 for k in range(n+m))
	for k,item in enumerate(case1):
		invertedindex[item[1][0]] = k
	for k,item in enumerate(case12):
		invertedindex[item[1][0]] = k+n1
		invertedindex[item[1][1]] = k+n1
	for k,item in enumerate(case2):
		invertedindex[item[1][0]] = k+n1+n2
@}
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
\subsubsection{Re-indexing of d-cells}
%-------------------------------------------------------------------------------

\paragraph{Return the single reordered pointset and the two $d$-cell arrays}
We are now finally ready to return two reordered LAR models defined over the same set \texttt{V} of vertices, and where (a) the vertex array \texttt{V} can be written as the union of three disjoint sets of points $C_1,C_{12},C_2$; (b) the $d$-cell array \texttt{CV1} is indexed over $C_1\cup C_{12}$; (b) the $d$-cell array \texttt{CV2} is indexed over $C_{12}\cup C_{2}$. 

The \texttt{vertexSieve} function will return the new reordered vertex set $V = (V_1 \cup V_2) \setminus (V_1 \cap V_2)$, the two renumbered $s$-cell sets \texttt{CV1} and \texttt{CV2}, and the size \texttt{len(case12)} of $V_1 \cap V_2$.
%-------------------------------------------------------------------------------
@D Return the single reordered pointset and the two $d$-cell arrays
@{
	V = [eval(p[0]) for p in case1] + [eval(p[0]) for p in case12] + [eval(
				p[0]) for p in case2]
	CV1 = [sorted([invertedindex[v] for v in cell]) for cell in CV1]
	CV2 = [sorted([invertedindex[v] for v in cell]) for cell in CV2]
	return V, CV1, CV2, len(case12)
@}
%-------------------------------------------------------------------------------


\subsubsection{Example of input with some coincident vertices}
In this example we give two very simple LAR representations of 2D cell complexes, with some coincident vertices, and go ahead to re-index the vertices, according to the method implemented by the function \texttt{vertexSieve}.

%-------------------------------------------------------------------------------
@o test/py/boolean2/test02.py
@{@< Initial import of modules @>
@< Import the module @(boolean2@) @>
V1 = [[1,1],[3,3],[3,1],[2,3],[2,1],[1,3]]
V2 = [[1,1],[1,3],[2,3],[2,2],[3,2],[0,1],[0,0],[2,0],[3,0]]
CV1 = [[0,3,4,5],[1,2,3,4]]
CV2 = [[3,4,7,8],[0,1,2,3,5,6,7]]
model1 = V1,CV1; model2 = V2,CV2
VIEW(STRUCT([ 
	COLOR(CYAN)(SKEL_1(STRUCT(MKPOLS(model1)))), 
	COLOR(RED)(SKEL_1(STRUCT(MKPOLS(model2)))) ]))
V, n1,n2,n12, B1,B2 = boolOps(model1,model2)
# VIEW(SKEL_1(STRUCT(MKPOLS((V, CV_un[:n1]+CV_int )))))
# VIEW(SKEL_1(STRUCT(MKPOLS((V, CV_un[n1-n12:]+CV_int )))))
@}
%-------------------------------------------------------------------------------

\paragraph{Example discussion} 
The aim of the \texttt{vertexSieve} function is twofold: (a) eliminate vertex duplicates before entering the main part of the Boolean algorithm; (b) reorder the input representations so that it becomes less expensive to check whether a 0-cell can be shared by both the arguments of a Boolean expression, so that its coboundaries must be eventually split. Remind that for any set it is:
\[
|A\cup B| = |A|+|B|-|A\cap B|.
\]
Let us notice that in the previous example
\[
|V| = |V_1 \cup V_2| = 12 \leq |V_1|+|V_2| = 6+9 = 15,
\]
and that 
\[
|V_1|+|V_2| - |V_1 \cup V_2| = 15 - 12 = 3 = |C_{12}| = |V_1 \cap V_2|,
\]
where $C_{12}$ is the subset of vertices with duplicated instances.
%-------------------------------------------------------------------------------
@D Output from \texttt{test/py/boolean2/test02.py}
@{V   = [[3.0,1.0],[2.0,1.0],[3.0,3.0],[1.0,1.0],[1.0,3.0],[2.0,3.0],
		 [3.0,2.0],[2.0,0.0],[2.0,2.0],[0.0,0.0],[3.0,0.0],[0.0,1.0]]
CV1 = [[3,5,1,4],[2,0,5,1]]
CV2 = [[8,6,7,10],[3,4,5,8,11,9,7]]
@}
%-------------------------------------------------------------------------------
Notice also that \texttt{V} has been reordered in three consecutive subsets $C_{1},C_{12},C_{2}$ such that \texttt{CV1} is indexed within $C_{1}\cup C_{12}$, whereas \texttt{CV2} is indexed within $C_{12}\cup C_{2}$. In our example we have  $C_{12}=\{\texttt{3,4,5}\}$: 
%-------------------------------------------------------------------------------
@D Reordering of vertex indexing of cells
@{
>>> sorted(CAT(CV1))
[0, 1, 1, 2, 3, 4, 5, 5]
>>> sorted(CAT(CV2))
[3, 4, 5, 6, 7, 7, 8, 8, 9, 10, 11]
@}
%-------------------------------------------------------------------------------
\paragraph{Cost analysis} 
Of course, this reordering after elimination of duplicate vertices will allow to perform a cheap $O(n)$ discovering of (Delaunay) cells whose vertices belong both to \texttt{V1} \emph{and} to \texttt{V2}. 
Actually, the \emph{same test} can be now used both when the vertices of the input arguments are all different, \emph{and} when they have some coincident vertices.
The total cost of such pre-processing, executed using dictionaries, is $O(n\ln n)$.


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
\section{Extracting divisor $d$-cells}
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

It is well-known that, in order to compute any Boolean operation with cell-decomposed arguments, either using a decompositive or a boundary scheme, the intersection of boundaries must be computed~\cite{Paoluzzi:1989:BAO:70248.70249}. In this section we develop the preparatory work for such a task, aiming to compute in the end the coboundary of boundary of each Boolean argument. We call the $d$-cells in such two sets, namely
\[
\Delta^1 = (\delta_{d-1}\circ\partial_d)(\Lambda^1_d)
\quad\mbox{and}\quad
\Delta^2 = (\delta_{d-1}\circ\partial_d)(\Lambda^2_d),
\]
as the \emph{divisors} of $\Lambda^2$ and $\Lambda^1$, respectively.


\paragraph{Coboundary of boundary}

The function \texttt{coboundaryOfBoundaryCells} returns the $d$-cells of the coboundary of boundary, providing as input the \texttt{BRC} representation of the $d$-cells and the $(d-1)$-cells. In other words, it computes the set of $d$-cells
\[
\Delta = (\delta_{d-1}\circ\partial_d)(\Lambda_d)
\]
where $\texttt{cells} := \texttt{BRC}(\Lambda_d)$, and $\texttt{facets} := \texttt{BRC}(\Lambda_{d-1})$.

%-------------------------------------------------------------------------------
@D Compute the $d$-cells in the coboundary of boundary
@{def coboundaryOfBoundaryCells(cells,facets):
    csrBoundaryMat = boundary(cells,facets)
    csrChain = totalChain(cells)
    csrBoundaryChain = matrixProduct(csrBoundaryMat, csrChain)
    for k,value in enumerate(csrBoundaryChain.data):
        if value % 2 == 0: csrBoundaryChain.data[k] = 0
    boundaryCells = [k for k,val in enumerate(csrBoundaryChain.data.tolist())
                               if val == 1]
    csrCoboundaryBoundaryChain = matrixProduct(csrBoundaryMat.T, csrBoundaryChain)
    coboundaryBoundaryCells = [k for k in csrCoboundaryBoundaryChain.indices.tolist()]
    return coboundaryBoundaryCells
@}
%-------------------------------------------------------------------------------


\paragraph{Boundary divisors}
The actual computation of the two sets of boundary divisors $\Delta^1$ and $\Delta^2$ is performed by the \texttt{boundarySuperCells} function below.

%-------------------------------------------------------------------------------
@D Compute the boundary divisors of both arguments
@{""" Second stage of Boolean operations """
def boundarySuperCells( V, CV ):
	FV = larSimplexFacets(CV)
	BSupCells = coboundaryOfBoundaryCells(CV,FV)
	return BSupCells
@}
%-------------------------------------------------------------------------------


By definition, a \emph{divisor} is a single $d$-cell and, in the present implementation, it is a $d$-simplex. \texttt{BSupCells} is the array of $d$-cells in the coboundary of boundary. They become \emph{divisors} of the (other argument's boundary) iff NONE $d$-star of their vertices does contain the $d$-cell itself.

%-------------------------------------------------------------------------------
@D Minimal covering of divisors
@{""" Minimal covering chains of divisors """
def minimalCovers(V,divisors,CV1,CV,setCV):
	covers = [selectIncidentChain( V, CV )(v) for cell in divisors for v in CV1[cell] ]
	print "\n covers =",covers
	return covers
@}
%-------------------------------------------------------------------------------




%-------------------------------------------------------------------------------
\subsection{Boundary computation}
%-------------------------------------------------------------------------------
The matrices of the boundary operators of the boolean arguments $\Lambda^1$ and $\Lambda^2$ are computed here as supported by the novel vertex set $V := V_1 \cup V_2$. Both the characteristic matrices $M_d$ and $M_{d-1}$ are needed to compute a $[\partial_d]$ matrix (see Reference~\cite{Dicarlo:2014:TNL:2543138.2543294}). Hence we start this section by computing the new basis of $(d-1)$-faces $\texttt{FV} := \texttt{CSR}(M_{d-1})$, and then compute the two subsets $B^1,B^2 \subset V$ of boundary vertices (upon the joint Delaunay complex \texttt{V}), where
\[
B^1 = [\mathcal{VF}^1]\, [\partial_d^1]\,\mathbf{1},
\quad\mbox{and}\quad
B^2 = [\mathcal{VF}^2]\, [\partial_d^2]\,\mathbf{1}.
\]
where $[\mathcal{VF}^1]^\top = \texttt{CSR(FV1)}$ and $[\mathcal{VF}^2]^\top = \texttt{CSR(FV2)}$,
and where \texttt{FV1} and \texttt{FV2} are the relations \emph{face-vertices} computed from the relation \texttt{CV} supported by the \emph{joint} Delaunay vertex set. 

\paragraph{Compute the boundary vertices of both arguments}
The two bases of $d$-cells, given as input to the \texttt{boundaryVertices} function below, were already renumbered. In other words, their vertices currently belong to the common Delaunay complex. Therefore, the subsequent calls to \texttt{larSimplexFacets} also return two sets of boundary facets, denoted as \texttt{BF1} and \texttt{BF2}, are supported by the Delaunay complex. \texttt{BV1} and \texttt{BV2} contain the boundary vertex indices of the input cells \texttt{CV1} and \texttt{CV2}.

%-------------------------------------------------------------------------------
@D Compute boundary vertices of both arguments
@{""" Second stage of Boolean operations """
def boundaryVertices( V, CV1,CV2 ):
	FV1 = larSimplexFacets(CV1)
	FV2 = larSimplexFacets(CV2)
	BF1 = boundaryCells(CV1,FV1)
	BF2 = boundaryCells(CV2,FV2)
	BV1 = list(set([ v for f in BF1 for v in FV1[f] ]))
	BV2 = list(set([ v for f in BF2 for v in FV2[f] ]))
	VIEW(STRUCT([ 
		COLOR(GREEN)(STRUCT(AA(MK)([V[v] for v in BV1]))), 
		COLOR(MAGENTA)(STRUCT(AA(MK)([V[v] for v in BV2]))) ]))
	return BV1, BV2
@}
%-------------------------------------------------------------------------------


\begin{figure}[htbp] %  figure placement: here, top, bottom, or page
   \includegraphics[height=0.245\linewidth,width=0.245\linewidth]{images/boolchain0a} 
   \includegraphics[height=0.245\linewidth,width=0.245\linewidth]{images/boolchain0b} 
   \includegraphics[height=0.245\linewidth,width=0.245\linewidth]{images/boolchain0c} 
   \includegraphics[height=0.245\linewidth,width=0.245\textwidth]{images/boolchain0d} 

   \includegraphics[height=0.245\linewidth,width=0.245\linewidth]{images/boolchain1a} 
   \includegraphics[height=0.245\linewidth,width=0.245\linewidth]{images/boolchain1b} 
   \includegraphics[height=0.245\linewidth,width=0.245\linewidth]{images/boolchain1c} 
   \includegraphics[height=0.245\linewidth,width=0.245\textwidth]{images/boolchain1d} 
   \caption{The selection of chains candidate to splitting, in both unstructured and structured meshes.}
   \label{fig:example}
\end{figure}


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
\section{Splitting argument chains}
\label{sec:splitting}
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

The \texttt{boolOps} function is the main procedure of the Boolean algorithm.
Its steps give an outline of the computations to be performed in sequence. The input LAR representations are first decomposed in their geometric and topological components, and embedded in the same space, so obtaining a global renumbering as \texttt{V}, \texttt{CV1}, \texttt{CV2}. Then the global Delaunay's triangulation \texttt{CV} and its set representation \texttt{setCV} are computed. In due order, the following steps are:  
\begin{enumerate}
\item 
the extraction of $\gamma^1,\gamma^2\subseteq\Lambda_d$, the coboundary of the boundaries of the Boolean arguments:
\[
\gamma^1 = (\delta_{d-1}\circ\partial_d) (\Lambda_d^1), \qquad
\gamma^2 = (\delta_{d-1}\circ\partial_d) (\Lambda_d^2);
\]  
\item 
the selection of the two invariant sub-chains, i.e.~of the subsets of $\Lambda_d^1, \Lambda_d^2$ that were not changed by the re-triangulation, and hence continue to contain some facets of the original boundaries:
\[
\gamma^{1'} \subseteq \gamma^1, \qquad
\gamma^{2'} \subseteq \gamma^2
\]
\end{enumerate}

\paragraph{The top-level algorithm}

%-------------------------------------------------------------------------------
@D Boolean subdivided complex
@{""" High level Boolean Application Programming Interface """
@< Minimal covering chains of divisors @>
def boolOps(lar1,lar2):
	(V1,CV1),(V2,CV2) = lar1,lar2
	n1,n2 = len(V1),len(V2)
	V, CV1, CV2, n12 = vertexSieve(lar1, lar2)
	CV = Delaunay(array(V)).vertices
	setCV = set([tuple(sorted(cell)) for cell in CV])
	@< Extraction of coboundary of boundary chains @>
	@< Computation of invariant supercells @>
	@< Invariant subsets of supercells @>
	@< Minimal covering of divisors @>
	@< Covers of divisors computation @>
	return V,n1,n2,n12, B1,B2
@}
%-------------------------------------------------------------------------------

\paragraph{Extraction of coboundary of boundary chains}

After the embedding of the two topologies in the same coordinate space and the Delaunay's re-triangulation, the two $\gamma^1,\gamma^2$ $d$-chains are computed by the following macro, and stored in the \texttt{BSupCells1} and \texttt{BSupCells2} variables.

%-------------------------------------------------------------------------------
@D Extraction of coboundary of boundary chains
@{# Second stage of Boolean algorithm
B1,B2 = boundaryVertices( V, CV1, CV2 )
# Extraction of coboundary of boundary chains
BSupCells1 = boundarySuperCells( V, CV1 )
BSupCells2 = boundarySuperCells( V, CV2 )

VIEW(STRUCT([ 
	COLOR(GREEN)(STRUCT(MKPOLS(((V,[CV1[c] for c in BSupCells1]))))), 
	COLOR(MAGENTA)(STRUCT(MKPOLS(((V,[CV2[c] for c in BSupCells2]))))) 
])) 
@}
%-------------------------------------------------------------------------------

\paragraph{Invariant subsets of supercells}

Then the sub-chains $\gamma^{1'},\gamma^{2'}$ (those not changed by the Boolean merging) are computed by the \texttt{invariantSuperCells} function, and stored in \texttt{invSupCells1} and \texttt{invSupCells2} variables, in order to compute, by chain difference, the \texttt{divisor1} and \texttt{divisor2} $d$-chains.

%-------------------------------------------------------------------------------
@D Computation of invariant supercells
@{def invariantSuperCells(V,BSupCells1,CV1,setCV):
	out = []
	for cell in BSupCells1:
		if tuple(CV1[cell]) in setCV: out += [cell]
	return out
@}
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
@D Invariant subsets of supercells
@{# Invariant subsets of supercells
invSupCells1 = invariantSuperCells(V,BSupCells1,CV1,setCV)
invSupCells2 = invariantSuperCells(V,BSupCells2,CV2,setCV)

divisor1 = set(BSupCells1).difference(invSupCells1)
divisor2 = set(BSupCells2).difference(invSupCells2)
VIEW(STRUCT([ 
	COLOR(GREEN)(STRUCT(MKPOLS(((V,[CV1[c] for c in divisor1]))))), 
	COLOR(MAGENTA)(STRUCT(MKPOLS(((V,[CV2[c] for c in divisor2]))))) 
]))
@}
%-------------------------------------------------------------------------------
\paragraph{Covers of divisors computation}
%-------------------------------------------------------------------------------
@D Covers of divisors computation
@{
# Covers of divisors computation
covers1 = minimalCovers(V,divisor1,CV1,CV,setCV)
covers2 = minimalCovers(V,divisor2,CV2,CV,setCV)	
VIEW(STRUCT([ 
	EXPLODE(1.2,1.2,1)(MKPOLS(((V,CV)))), 
	COLOR(GREEN)(EXPLODE(1.2,1.2,1)(MKPOLS(((V,[CV[c] for c in CAT(covers1)]))))), 
	COLOR(MAGENTA)(EXPLODE(1.2,1.2,1)(MKPOLS(((V,[CV[c] for c in CAT(covers2)])))))
	 ]))
@}
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
\subsection{Matching cells in $\Sigma_\cap$ with spanning chains in $\Lambda^1$, $\Lambda^2$}
%-------------------------------------------------------------------------------

The next step of the \emph{Boolean Chains} algorithm retrieves and compare the chains of $d$-cells incident on $\texttt{V}(\partial_d(\Lambda^1)), \texttt{V}(\partial_d(\Lambda^2))$ and on $\texttt{V}(\partial_d(\Sigma^\cap))$. A filtering step working in $O(n)$ is performed in regard to this.

\paragraph{Filtering cells incident on a subset of vertices}

The selection of the $d$-chain incident on a $0$-chain (subset of vertices) can be executed either in time roughly proportional to the output, via a SpMV multiplication of the $\texttt{CV} := \texttt{CSR}(M_d)$ matrix times the (binary) coordinate representation of the $0$-chain, or in time proportional to the input, via a more traditional select filtering of cells in a \texttt{CV} table, with respect to the subset of vertices.


%-------------------------------------------------------------------------------
@D Select cells incident on vertices
@{""" Select the $d$-chain incident on a $0$-chain """
def selectIncidentChain( V, CV ):
	def selectIncidentChain0( v ):
		return [k for k in range(len(CV)) if v in CV[k] ]
	return selectIncidentChain0
@}
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
\subsection{Splitting cells}
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
\subsection{Keeping cell dictionaries updated}
%-------------------------------------------------------------------------------
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
\section{Boolean outputs computations}
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
\section{Export the boolean module}
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

The \texttt{boolean.py} module is exported to the library \texttt{lar-cc/lib}. Therefore many of the macros developed in this module are expanded and written to an external file.

%------------------------------------------------------------------
@o lib/py/boolean2.py
@{""" Module with Boolean operators using chains and CSR matrices """
@< Initial import of modules @>
@< Symbolic utility to represent points as strings @>
@< Affine transformations of $d$-points @>
@< Generation of $n$ random points in the unit $d$-disk @>
@< Generation of $n$ random points in the standard $d$-cuboid @>
@< Triangulation of random points @>
@< Boolean subdivided complex @>
@< High-level boolean operations @>
@< Place the vertices of Boolean arguments in a common space @>
@< Show vertices of arguments @>
@< Compute boundary vertices of both arguments @>
@< Select cells incident on vertices @>
@< Compute the $d$-cells in the coboundary of boundary @>
@< Compute the boundary divisors of both arguments @>
@< Visualization of a subset of cells @>
@}
%------------------------------------------------------------------

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
\section{Tests}
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%-------------------------------------------------------------------------------
\subsection{Generation of random data}

We found useful to drive the development of new modules using randomly generated data, so that every upcoming execution of the developed algorithms is naturally driven to be challenged by different data.

\subsubsection{Testing the main algorithm}

\paragraph{Write the test executable file}




\subsubsection{Lowest-level space generation procedures}

\paragraph{Random points in unit disk} 
First we generate a  set of $n$ random points in the unit $D^d$ disk centred on the origin, to be subsequently used to generate a random Delaunay complex of variable granularity.

%------------------------------------------------------------------
@D Generation of $n$ random points in the unit $d$-disk 
@{def randomPointsInUnitCircle(n=200,d=2, r=1):
	points = random.random((n,d)) * ([2*math.pi]+[1]*(d-1))
	return [[SQRT(p[1])*COS(p[0]),SQRT(p[1])*SIN(p[0])] for p in points]
	## TODO: correct for $d$-sphere

if __name__=="__main__":
	VIEW(STRUCT(AA(MK)(randomPointsInUnitCircle()))) 
@}
%------------------------------------------------------------------

\paragraph{Random points in the standard $d$-cuboid} 
A set of $n$ random $d$-points is then generated within the standard $d$-cuboid, i.e.~withing the $d$-dimensional interval with a vertex on the origin.

%------------------------------------------------------------------
@D Generation of $n$ random points in the standard $d$-cuboid 
@{def randomPointsInUnitCuboid(n=200,d=2):
	return random.random((n,d)).tolist()

if __name__=="__main__":
	VIEW(STRUCT(AA(MK)(randomPointsInUnitCuboid()))) 
@}
%------------------------------------------------------------------



\paragraph{Triangulation of random points} The Delaunay triangulation of \texttt{randomPointsInUnitCircle} is generated by the following macro.


%------------------------------------------------------------------
@D Triangulation of random points
@{from scipy.spatial import Delaunay
def randomTriangulation(n=200,d=2,out='disk'):
	if out == 'disk':
		V = randomPointsInUnitCircle(n,d)
	elif out == 'cuboid':
		V = randomPointsInUnitCuboid(n,d)
	CV = Delaunay(array(V)).vertices
	model = V,CV
	return model

if __name__=="__main__":
	from lar2psm import *
	VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(model)))
@}
%------------------------------------------------------------------

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
\subsection{Disk saving of test data}
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
\subsection{Algorithm execution}
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
\subsection{Unit tests}
%-------------------------------------------------------------------------------

\subsubsection{First Boolean stage}

Some unit tests of the first Boolean stage are discussed in the following. They are mainly aimed to check a correct execution of the filtering of common vertices with renumbering of the union set of vertices, and to the consequential redefinition of the $d$-cell basis.

\paragraph{Union of 2D non-structured grids}

%------------------------------------------------------------------
@o test/py/boolean2/test01.py
@{""" Union of 2D non-structured grids """
@< Initial import of modules @>
@< Import the module @(boolean2@) @>
@< Import the module @(lar2psm@) @>
@< Import the module @(myfont@) @>
model1 = randomTriangulation(1000,2,'disk')
V1,CV1 = model1
VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(model1)+cellNames(model1,CV1,MAGENTA)))
model2 = randomTriangulation(1000,2,'cuboid')
V2,CV2 = model2
V2 = larScale( [2,2])(V2)
model2 = V2,CV2 
VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(model2)+cellNames(model2,CV2,RED)))
V, n1,n2,n12, B1,B2 = boolOps(model1,model2)
@< Visualization of first Boolean step @>
@}
%------------------------------------------------------------------

\paragraph{Union of 3D non-structured grids}

%------------------------------------------------------------------
@o test/py/boolean2/test05.py
@{""" Union of 3D non-structured grids """
@< Initial import of modules @>
@< Import the module @(boolean2@) @>
@< Import the module @(lar2psm@) @>
@< Import the module @(myfont@) @>
model1 = randomTriangulation(100,3,'cuboid')
V1,CV1 = model1
V1 = larScale( [2,2,2])(V1)
V1 = larTranslate( [-1,-1,-1])(V1)
model1 = V1,CV1 
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model1)+cellNames(model1,CV1,MAGENTA)))
model2 = randomTriangulation(100,3,'cuboid')
V2,CV2 = model2
V2 = larScale( [2,2,2])(V2)
model2 = V2,CV2 
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model2)+cellNames(model2,CV2,RED)))
V, n1,n2,n12, B1,B2 = boolOps(model1,model2)
@< Visualization of first Boolean step @>
@}
%------------------------------------------------------------------


\paragraph{Union of structured grids}

%------------------------------------------------------------------
@o test/py/boolean2/test04.py
@{""" test program for the boolean module """
@< Initial import of modules @>
@< Import the module @(boolean2@) @>
blue = larSimplexGrid([30,60])
V2,CV2 = larSimplexGrid([70,40])
V2 = larTranslate( [.5,.5])(V2)
red = V2,CV2
VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(blue) ))
VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(red) ))
V, CV1, CV2, n12 = vertexSieve(red,blue)
V, n1,n2,n12, B1,B2 = boolOps(red,blue)
CV = Delaunay(array(V)).vertices
@< Visualization of first Boolean step @>
@}
%------------------------------------------------------------------


model = checkModel(larHollowCyl(0.8,1,1,angle=PI/4)([12,2,2]))
VIEW(STRUCT(MKPOLS(model)))
model = checkModel(larHollowSphere(0.8,1,PI/6,PI/4)([6,12,2]))


\paragraph{Union of structured grids}

%------------------------------------------------------------------
@o test/py/boolean2/test06.py
@{""" test program for the boolean module """
@< Initial import of modules @>
from boolean2 import *
from mapper import *
blue = larHollowCyl(0.8,1,1,angle=PI/4)([6,2,5])
red = larHollowSphere(0.8,1,PI/6,PI/4)([6,12,2])
VIEW(EXPLODE(1.5,1.5,1)(MKPOLS(blue) + MKPOLS(red) ))
V, CV1, CV2, n12 = vertexSieve(red,blue)
V, n1,n2,n12, B1,B2 = boolOps(red,blue)
CV = Delaunay(array(V)).vertices
@< Visualization of first Boolean step @>
@}
%------------------------------------------------------------------





\paragraph{Boolean operations with general LAR cells}

%------------------------------------------------------------------
@o test/py/boolean2/test03.py
@{""" test example with general LAR cells for the boolean module """
@< Initial import of modules @>
@< Import the module @(boolean2@) @>
V1 = [[0,5],[6,5],[0,2],[3,2],[6,2],[0,0],[3,0],[3,-2],[6,-2]]
CV1 = [[2,3,5,6],[0,1,2,3,4],[3,4,6,7,8]]
blue = V1,CV1
V2 = [[3,6],[7,6],[0,5],[3,5],[3,4],[7,4],[3,2],[7,2],[0,0],[3,0],[6,0],[6,2]]
CV2 = [[0,1,3,4,5],[2,3,4,6,8,9],[6,9,10,11],[4,5,6,7,11]]
red = V2,CV2
V, n1,n2,n12, B1,B2 = boolOps(red,blue)
CV = Delaunay(array(V)).vertices
@< Visualization of first Boolean step @>
@}
%------------------------------------------------------------------

\subsection{Visualization of the Boolean algorithm}

\paragraph{Display of colored cell numbers}

%------------------------------------------------------------------
@D Visualization of a subset of cells
@{def cellNames(model,cells, color=BLACK):
	V,CV= model
	print "\n CV =",CV
	print "\n cells =",cells
	texts = []
	for k,cell in enumerate(cells):
		centroid = CCOMB([V[v] for v in cell])
		print "centroid =",centroid
		d = len(centroid)
		texts += [T(range(1,d+1))(centroid)(S(range(1,d+1))([0.005 for k in range(d)])(TEXT(str(k))))]
	return AA(COLOR(color))(texts)
@}
%------------------------------------------------------------------


\paragraph{First step of visualization}

%------------------------------------------------------------------
@D Visualization of first Boolean step
@{""" Visualization of first Boolean step  """
if n12==0:
	hpc0 = STRUCT([ COLOR(RED)(EXPLODE(1.5,1.5,1)(AA(MK)(V[:n1-n12]) )), 
				COLOR(CYAN)(EXPLODE(1.5,1.5,1)(AA(MK)(V[n1:]) )) ])
else:
	hpc0 = STRUCT([ COLOR(RED)(EXPLODE(1.5,1.5,1)(AA(MK)(V[:n1-n12]) )), 
				COLOR(CYAN)(EXPLODE(1.5,1.5,1)(AA(MK)(V[n1:]) )), 
				COLOR(WHITE)(EXPLODE(1.5,1.5,1)(AA(MK)(V[n1-n12:n1]) )) ])

# hpc1 = COLOR(RED)(EXPLODE(1.5,1.5,1)(MKPOLS((V,CV_un)) ))
# hpc2 = COLOR(CYAN)(EXPLODE(1.5,1.5,1)(MKPOLS((V,CV_int)) ))
# VIEW(STRUCT([hpc0, hpc1, hpc2]))
@}
%------------------------------------------------------------------

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
\appendix
\section{Utility functions}

%------------------------------------------------------------------
@D Initial import of modules
@{from pyplasm import *
from scipy import *
import os,sys

""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
@< Import the module @(lar2psm@) @>
@< Import the module @(simplexn@) @>
@< Import the module @(larcc@) @>
@< Import the module @(largrid@) @>
@< Import the module @(myfont@) @>
@< Import the module @(mapper@) @>
@}
%------------------------------------------------------------------

%------------------------------------------------------------------
@D Import the module
@{import @1
from @1 import *
@}
%------------------------------------------------------------------

\paragraph{Affine transformations of points} Some primitive maps of points to points are given in the following, including translation, rotation and scaling of array of points via direct transformation of their coordinates.

%------------------------------------------------------------------
@D Affine transformations of $d$-points
@{def translatePoints (points, tvect):
	return [VECTSUM([p,tvect]) for p in points]

def rotatePoints (points, angle):		# 2-dimensional !! TODO: n-dim
	a = angle
	return [[x*COS(a)-y*SIN(a), x*SIN(a)+y*COS(a)] for x,y in points]

def scalePoints (points, svect):
	return [AA(PROD)(TRANS([p,svect])) for p in points]
@}
%------------------------------------------------------------------

\subsection{Numeric utilities}

A small set of utilityy functions is used to transform a point representation as array of coordinates into a string of fixed format to be used as point key into python dictionaries.

%------------------------------------------------------------------
@D Symbolic utility to represent points as strings
@{""" TODO: use package Decimal (http://docs.python.org/2/library/decimal.html) """
ROUND_ZERO = 1E-07
def round_or_zero (x,prec=7):
	"""
	Decision procedure to approximate a small number to zero.
	Return either the input number or zero.
	"""
	def myround(x):
		return eval(('%.'+str(prec)+'f') % round(x,prec))
	xx = myround(x)
	if abs(xx) < ROUND_ZERO: return 0.0
	else: return xx

def prepKey (args): return "["+", ".join(args)+"]"

def fixedPrec(value):
	if abs(value - int(value))<ROUND_ZERO: value = int(value)
	out = ('%0.7f'% value).rstrip('0')
	if out == '-0.': out = '0.'
	return out
	
def vcode (vect): 
	"""
	To generate a string representation of a number array.
	Used to generate the vertex keys in PointSet dictionary, and other similar operations.
	"""
	return prepKey(AA(fixedPrec)(vect))
@}
%------------------------------------------------------------------




\bibliographystyle{amsalpha}
\bibliography{boolean2}

\end{document}
