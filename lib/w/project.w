\documentclass[11pt,oneside]{article}	%use"amsart"insteadof"article"forAMSLaTeXformat
\usepackage{geometry}		%Seegeometry.pdftolearnthelayoutoptions.Therearelots.
\geometry{letterpaper}		%...ora4paperora5paperor...
%\geometry{landscape}		%Activateforforrotatedpagegeometry
%\usepackage[parfill]{parskip}		%Activatetobeginparagraphswithanemptylineratherthananindent
\usepackage{graphicx}				%Usepdf,png,jpg,orepsßwithpdflatex;useepsinDVImode
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

\title{Package per il mio progetto
\footnote{This document is part of the \emph{Linear Algebraic Representation with CoChains} (LAR-CC) framework~\cite{cclar-proj:2013:00}. \today}
}
\author{Nome1 Cognome1 \and Nome2 Cognome2}
%\date{}							%Activatetodisplayagivendateornodate

\begin{document}
\maketitle
\nonstopmode

\begin{abstract}
aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa aaaa aaa aaaa 
\end{abstract}

\tableofcontents

\section{La mia prima sezione}

\subsection{La mia prima sottosezionesezione}

\subsection{La mia seconda sottosezionesezione}


\section{La mia seconda sezione}


\section{La mia terza sezione}

\bibliographystyle{amsalpha}
\bibliography{project}

\end{document}
