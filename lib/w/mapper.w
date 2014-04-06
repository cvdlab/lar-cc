\documentclass[11pt,oneside]{article}	%use"amsart"insteadof"article"forAMSLaTeXformat
\usepackage{geometry}		%Seegeometry.pdftolearnthelayoutoptions.Therearelots.
\geometry{letterpaper}		%...ora4paperora5paperor...
%\geometry{landscape}		%Activateforforrotatedpagegeometry
%\usepackage[parfill]{parskip}		%Activatetobeginparagraphswithanemptylineratherthananindent
\usepackage{graphicx}				%Usepdf,png,jpg,orepsßwithpdflatex;useepsinDVImode
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

\title{Domain mapping with LAR
\footnote{This document is part of the \emph{Linear Algebraic Representation with CoChains} (LAR-CC) framework~\cite{cclar-proj:2013:00}. \today}
}
\author{Alberto Paoluzzi}
%\date{}							%Activatetodisplayagivendateornodate

\begin{document}
\maketitle
\nonstopmode

\begin{abstract}
In this module a first implementation (no optimisations) is done of the \texttt{LARMAP} operator, reproducing the behaviour of the plasm \texttt{MAP} primitive, but with better handling of the topology, including the sewing of decomposed (simplicial domains) about their possible sewing.
\end{abstract}

\tableofcontents

%===============================================================================
\section{Domain decomposition}
%===============================================================================

\paragraph{Standard and scaled decomposition of unit domain}

%-------------------------------------------------------------------------------
@d Generate a simplicial decomposition ot the $[0,1]^d$ domain
@{def larDomain(shape):
	V,CV = larSimplexGrid(shape)
	V = scalePoints(V, [1./d for d in shape])
	return V,CV

def larIntervals(shape):
	def larIntervals0(size):
		V,CV = larDomain(shape)
		V = scalePoints(V, [scaleFactor for scaleFactor in size])
		return V,CV
	return larIntervals0
@}
%-------------------------------------------------------------------------------

%===============================================================================
\section{Embedding via coordinate functions}
%===============================================================================

%-------------------------------------------------------------------------------
@D Primitive mapping function 
@{def larMap(coordFuncs):
	def larMap0(domain):
		V,CV = domain
		V = TRANS(CONS(coordFuncs)(V))
		return V,CV
	return larMap0
@}
%-------------------------------------------------------------------------------

%===============================================================================
\section{Primitive objets}
%===============================================================================
%-------------------------------------------------------------------------------
\subsection{1D primitives}
%-------------------------------------------------------------------------------

\paragraph{Curved line}
%-------------------------------------------------------------------------------
@D aaaa
@{

@}
%-------------------------------------------------------------------------------

\paragraph{Circle}
%-------------------------------------------------------------------------------
@D Circumference of unit radius
@{def larCircle(radius=1.):
	def larCircle0(shape=36):
		domain = larIntervals([shape])([2*PI])
		V,CV = domain
		x = lambda coords : [radius*COS(p[0]) for p in V]
		y = lambda coords : [radius*SIN(p[0]) for p in V]
		mapping = [x,y]
		return larMap(mapping)(domain)
	return larCircle0

@}
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
\subsection{2D primitives}
%-------------------------------------------------------------------------------

\paragraph{Disk}
%-------------------------------------------------------------------------------
@D aaaa
@{

@}
%-------------------------------------------------------------------------------

\paragraph{Sphere surface}
%-------------------------------------------------------------------------------
@D aaaa
@{

@}
%-------------------------------------------------------------------------------

\paragraph{Cylinder surface}
%-------------------------------------------------------------------------------
@D aaaa
@{

@}
%-------------------------------------------------------------------------------

\paragraph{Torus surface}
%-------------------------------------------------------------------------------
@D aaaa
@{

@}
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
\subsection{3D primitives}
%-------------------------------------------------------------------------------

\paragraph{Ball}
%-------------------------------------------------------------------------------
@D aaaa
@{

@}
%-------------------------------------------------------------------------------

\paragraph{Solid cylinder}
%-------------------------------------------------------------------------------
@D aaaa
@{

@}
%-------------------------------------------------------------------------------

\paragraph{Solid torus}
%-------------------------------------------------------------------------------
@D aaaa
@{

@}
%-------------------------------------------------------------------------------

%===============================================================================
\section{Exporting the library}
%===============================================================================
%-------------------------------------------------------------------------------
@O lib/py/mapper.py
@{""" Mapping functions and primitive objects """
@< Initial import of modules @>
@< Generate a simplicial decomposition ot the $[0,1]^d$ domain @>
@< Primitive mapping function @>
@< Basic tests of mapper module @>
@< Circumference of unit radius @>
@}
%-------------------------------------------------------------------------------
%===============================================================================
\section{Examples}
%===============================================================================
%===============================================================================
\section{Tests}
%===============================================================================

	
%-------------------------------------------------------------------------------
@D Basic tests of mapper module
@{if __name__=="__main__":
	V,EV = larDomain([5])
	VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,EV))))
	V,EV = larIntervals([24])([2*PI])
	VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,EV))))
		
	V,FV = larDomain([5,3])
	VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))
	V,FV = larIntervals([36,3])([2*PI,1.])
	VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,FV))))
		
	V,CV = larDomain([5,3,1])
	VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))
	V,CV = larIntervals([36,2,3])([2*PI,1.,1.])
	VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))
@}
%-------------------------------------------------------------------------------

\paragraph{Circumference of unit radius}

%-------------------------------------------------------------------------------
@O test/py/mapper/test01.py
@{""" Circumference of unit radius """
@< Initial import of modules @>
import mapper
from mapper import larCircle
model = larCircle(1)()
VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
@}
%-------------------------------------------------------------------------------


%===============================================================================
\appendix
\section{Utility functions}
%===============================================================================

%-------------------------------------------------------------------------------
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
@< Import the module @(boolean2@) @>
@}
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
@D Import the module
@{import @1
from @1 import *
@}
%-------------------------------------------------------------------------------

\bibliographystyle{amsalpha}
\bibliography{mapper}

\end{document}
