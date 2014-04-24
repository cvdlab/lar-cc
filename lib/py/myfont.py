##%-----------------------------------------------------------------%
##%                                                                 %
##% Source codes extracted from:                                    %
##% A. Paoluzzi, Geometric Programming for Computer Aided Design,   %
##% J. Wiley & Sons Ltd, Chichester, UK, 2002                       %
##%                                                                 %
##%-----------------------------------------------------------------%

from pyplasm import *

##%----------------------------------------------------------------%
##%---Font as sequence of polyhedral complexes---------------------%
##%---Drawable ASCII subset [32-126]-------------------------------%
##%----------------------------------------------------------------%

ascii32 = Plasm.mkpol(2,[0.0,0.0],[[1]])
Xpol = COMP([STRUCT,AA(POLYLINE)])
ascii46 = Xpol([ [[2,0],[2,0.5],[1.5,0.5],[1.5,0],[2,0]] ])
ascii33 = STRUCT([ Xpol([ [[1.75,1.75],[1.75,5.5]] ]), T(2)(0.5)(ascii46) ])
ascii39 = STRUCT([ Xpol([ [[1,4],[2,5],[2,5.5]] ]), T(2)(5.5)(ascii46) ])
ascii34 = STRUCT([ ascii39, T(1)(1.0)(ascii39) ])
ascii61 = Xpol([ [[1,2.5],[3,2.5]],[[1,3.5],[3,3.5]] ])
ascii35 = STRUCT([ ascii61, Xpol([ [[1.25,1.75],[1.75,4]],[[2.25,1.75],[2.75,4]] ]) ])
ascii38 = Xpol([ [[4,1],[3,0],[1,0],[0,1],[0,2],[1,3],[2,3],[3,4],[3,5],[2,6],[1,5],[1,4],[3,2],[4,2]] ])
ascii40 = Xpol([ [[2,0],[1,1],[0.5,3],[1,5],[2,6]] ])
ascii41 = Xpol([ [[2,0],[3,1],[3.5,3],[3,5],[2,6]] ])
ascii43 = Xpol([ [[1,3],[3,3]],[[2,2],[2,4]] ])
ascii45 = Xpol([ [[1,3],[3,3]] ])
ascii47 = Xpol([ [[1,0],[3,6]] ])
ascii49 = Xpol([ [[0,4],[2,6],[2,0]],[[0,0],[4,0]] ])
ascii50 = Xpol([ [[0,4],[0,5],[1,6],[3,6],[4,5],[4,4],[0,0],[4,0]] ])
ascii51 = Xpol([ [[0,6],[4,6],[2,4],[4,2],[4,1],[3,0],[1,0],[0,1],[0,2]] ])
ascii52 = Xpol([ [[4,1],[0,1],[4,6],[4,0]] ])
ascii53 = Xpol([ [[4,6],[0,6],[0,4],[3,4],[4,3],[4,1],[3,0],[1,0],[0,1],[0,2]] ])
ascii54 = Xpol([ [[4,6],[1,6],[0,5],[0,1],[1,0],[3,0],[4,1],[4,3],[3,4],[1,4],[0,3]] ])
ascii55 = Xpol([ [[0,5],[0,6],[4,6],[0,0]] ] )
ascii57 = Xpol([ [[0,0],[3,0],[4,1],[4,5],[3,6],[1,6],[0,5],[0,3],[1,2],[3,2],[4,3]] ])
ascii60 = Xpol([ [[3,6],[0,3],[3,0]] ] )
ascii62 = Xpol([ [[1,6],[4,3],[1,0]] ] )
ascii64 = Xpol([ [[4,0],[1,0],[0,1],[0,3],[1,4],[3,4],[4,3],[4,1],[2,1],[1,2],[2,3],[3,2],[2,1]] ] )
ascii65 = Xpol([ [[0,0],[0,5],[1,6],[3,6],[4,5],[4,0]], [[0,2],[4,2]] ])
ascii66 = Xpol([ [[0,0],[0,6],[3,6],[4,5],[4,4],[3,3],[4,2], [4,1],[3,0],[0,0]],[[0,3],[3,3]]  ])
ascii67 = Xpol([ [[4,1],[3,0],[1,0],[0,1], [0,5],[1,6],[3,6],[4,5]] ])
ascii68 = Xpol([ [[0,0],[0,6],[3,6],[4,5],[4,1],[3,0],[0,0]] ])
ascii69 = Xpol([ [[4,0],[0,0],[0,6],[4,6]],[[0,3],[3,3]] ])
ascii70 = Xpol([ [[0,0],[0,6],[4,6]],[[0,3],[3,3]] ])
ascii71 = Xpol([ [[2,3],[4,3],[4,1],[3,0],[1,0],[0,1], [0,5],[1,6],[3,6],[4,5]] ])
ascii72 = Xpol([ [[0,0],[0,6]],[[4,6],[4,0]],[[0,3],[4,3]] ])
ascii73 = Xpol([ [[2,0],[2,6]],[[1,0],[3,0]],[[1,6],[3,6]] ])
ascii74 = Xpol([ [[0,1],[1,0],[2,0],[3,1],[3,6]],[[2,6],[4,6]] ])
ascii75 = Xpol([ [[4,6],[0,3],[4,0]],[[0,0],[0,6]] ])
ascii76 = Xpol([ [[4,0],[0,0],[0,6]] ])
ascii77 = Xpol([ [[0,0],[0,6],[2,4],[4,6],[4,0]] ])
ascii78 = Xpol([ [[0,0],[0,6],[4,2]],[[4,0],[4,6]] ])
ascii79 = Xpol([ [[4,1],[3,0],[1,0],[0,1], [0,5],[1,6],[3,6],[4,5],[4,1]] ])
ascii80 = Xpol([ [[0,0],[0,6],[3,6],[4,5],[4,3],[3,2],[0,2]] ])
ascii81 = Xpol([ [[4,1],[3,0],[1,0],[0,1], [0,5],[1,6],[3,6],[4,5],[4,1]],[[3,1],[4,0]] ])
ascii82 = Xpol([ [[0,0],[0,6],[3,6],[4,5],[4,3],[3,2],[0,2]], [[3,2],[4,0]] ])
ascii83 = Xpol([ [[0,1],[1,0],[3,0],[4,1],[4,2],[3,3],[1,3],[0,4],[0,5],[1,6],[3,6],[4,5]] ])
ascii84 = Xpol([ [[2,0],[2,6]],[[0,6],[4,6]] ])
ascii85 = Xpol([ [[0,6],[0,1],[1,0],[3,0],[4,1],[4,6]] ])
ascii86 = Xpol([ [[0,6],[2,0],[4,6]] ])
ascii87 = Xpol([ [[0,6],[0,3],[1,0],[2,3],[3,0],[4,3],[4,6]] ])
ascii88 = Xpol([ [[0,0],[4,6]],[[0,6],[4,0]] ])
ascii89 = Xpol([ [[0,6],[2,2],[4,6]],[[2,2],[2,0]] ])
ascii90 = Xpol([ [[0,6],[4,6],[0,0],[4,0]] ])
ascii91 = Xpol([ [[2,0],[1,0],[1,6],[2,6]] ])
ascii92 = Xpol([ [[1,6],[3,0]] ])
ascii93 = Xpol([ [[2,0],[3,0],[3,6],[2,6]] ])
ascii94 = Xpol([ [[1,5],[2,6],[3,5]] ])
ascii95 = Xpol([ [[1,0],[4,0]] ])
ascii99 = Xpol([ [[4,1],[3,0],[1,0],[0,1],[0,2],[1,3],[3,3],[4,2]] ])
ascii101 = Xpol([ [[4,0],[1,0],[0,1],[0,2],[1,3],[3,3],[4,2],[4,1],[0,1]] ])
ascii102 = Xpol([ [[4,3],[4,4],[3,5],[2,5],[1,4],[1,0]],[[0,1],[2,1]] ])
ascii104 = Xpol([ [[4,0],[4,2],[3,3],[1,3],[0,2]],[[0,0],[0,5],[1,5]] ])
ascii107 = Xpol([ [[0,0],[1,0],[1,3],[0,3]],[[4,0],[2,0],[1,1],[3,3],[4,3]] ])
ascii108 = Xpol([ [[2,0],[2,5],[1,5]],[[1,0],[3,0]] ])
ascii109 = Xpol([ [[4,0],[4,3],[2,2]],[[2,0],[2,3],[0,2]],[[0,0],[0,3]] ])
ascii110 = Xpol([ [[3,0],[3,3],[1,2]],[[1,0],[1,3]] ])
ascii111 = Xpol([ [[4,1],[3,0],[1,0],[0,1],[0,2],[1,3],[3,3],[4,2],[4,1]] ])
ascii114 = Xpol([ [[0,0],[2,0]],[[1,0],[1,3]],[[1,2],[2,3],[3,3],[4,2]] ])
ascii115 = Xpol([ [[0,0],[4,0],[3,1],[1,1],[0,2],[1,3],[3,3],[4,2]] ])
ascii116 = Xpol([ [[1,0],[3,0]],[[2,0],[2,5]],[[2,4],[3,4]] ])
ascii117 = Xpol([ [[0,3],[1,3],[1,1],[2,0],[3,0],[4,1],[4,3]] ])
ascii118 = Xpol([ [[0,3],[1,0],[3,3],[4,3]] ])
ascii119 = Xpol([ [[0,3],[0,2],[1,0],[2,2],[3,0],[4,2],[4,3]] ])
ascii120 = Xpol([ [[0,3],[1,3],[4,0]],[[1,0],[4,3]] ])
ascii121 = Xpol([ [[0,3],[1,3],[2.5,1.5]],[[0,0],[1,0],[4,3]] ])
ascii122 = Xpol([ [[0,2],[0,3],[3,3],[0,0],[3,0],[4,1]] ])
ascii123 = Xpol([ [[2.5,6.5],[2,6],[2,3.5],[1.5,3],[2,2.5],[2,0],[2.5,-0.5]] ])
ascii124 = Xpol([ [[2,0],[2,5]] ])
ascii125 = Xpol([ [[1.5,6.5],[2,6],[2,3.5],[2.5,3],[2,2.5],[2,0],[1.5,-0.5]] ])
ascii126 = Xpol([ [[1,5],[1.75,5.5],[2.75,5],[3.5,5.5]] ])
ascii36 = STRUCT([ ascii83, Xpol([ [[2,-0.5],[2,6.5]] ] )])
ascii37 = STRUCT([ T(1)(0.5)(ascii46), T(2)(5.5)(ascii46), Xpol([ [[1,0],[3,6]] ] )])
ascii39 = STRUCT([ Xpol([ [[1,4],[2,5],[2,5.5]] ]), T(2)(5.5)(ascii46) ])
ascii42 = STRUCT([ ascii43, COMP([ T([1,2])([2,3]), R([1,2])(PI/4), T([1,2])([-2,-3]) ])(ascii43) ])
ascii44 = T(2)(-5)(ascii39)
ascii48 = STRUCT([ ascii79, Xpol([ [[0,1],[4,5]] ]) ])
ascii56 = COMP([STRUCT, CONS([ID,T(2)(3)]), Xpol])([ [[1,0],[3,0],[4,1],[4,2],[3,3],[1,3],[0,2],[0,1],[1,0]] ])           
ascii58 = COMP([STRUCT, CONS([T(2)(1),T(2)(3)])])(ascii46) 
ascii59 = COMP([STRUCT, CONS([ COMP([T(2)(3), S1]), COMP([ T(2)(0.5), S2 ]) ])])([ascii46, ascii44]) 
ascii63 = STRUCT([ ascii46, Xpol([ [[1.75,1],[1.75,2.75],[3,4],[3,5],[2,6],[1,6],[0,5],[0,4]] ]) ])
ascii96 = STRUCT([ T([1,2])([0.5,4])(ascii46), Xpol([ [[2,4.5],[2,5],[3,6]] ]) ])
ascii97 = STRUCT([ ascii111, Xpol([ [[4,0],[4,3]] ] )])
ascii98 = STRUCT([ ascii111, Xpol([ [[0,0],[0,5],[1,5]] ] )])
ascii100 = STRUCT([ ascii111, Xpol([ [[4,0],[4,5],[3,5]] ] )])
ascii103 = STRUCT([ ascii111, Xpol([ [[4,1],[4,0],[3,-1],[1,-1],[0,0]] ] )])
ascii105 = STRUCT([ T([1,2])([0.25,3.75])(ascii46), Xpol([ [[1,0],[3,0]],[[1,3],[3,3]],[[2,0],[2,3]] ]) ])
ascii106 = STRUCT([ T([1,2])([0.25,3.75])(ascii46), Xpol([ [[1,3],[3,3]],[[2,3],[2,0],[1,-1],[0,0]] ]) ])
ascii112 = STRUCT([ ascii111, Xpol([ [[0,3],[0,-1]] ]) ])
ascii113 = STRUCT([ ascii111, Xpol([ [[4,3],[4,-1]] ]) ])

MyFont = [ 
	ascii32, ascii33, ascii34, ascii35, ascii36, ascii37, ascii38, ascii39, ascii40, ascii41, 
	ascii42, ascii43, ascii44, ascii45, ascii46, ascii47, ascii48, ascii49, ascii50, ascii51, 
	ascii52, ascii53, ascii54, ascii55, ascii56, ascii57, ascii58, ascii59, ascii60, ascii61, 
	ascii62, ascii63, ascii64, ascii65, ascii66, ascii67, ascii68, ascii69, ascii70, ascii71, 
	ascii72, ascii73, ascii74, ascii75, ascii76, ascii77, ascii78, ascii79, ascii80, ascii81, 
	ascii82, ascii83, ascii84, ascii85, ascii86, ascii87, ascii88, ascii89, ascii90, ascii91, 
	ascii92, ascii93, ascii94, ascii95, ascii96, ascii97, ascii98, ascii99, ascii100, ascii101, 
	ascii102, ascii103, ascii104, ascii105, ascii106, ascii107, ascii108, ascii109, ascii110, ascii111,
	ascii112, ascii113, ascii114, ascii115, ascii116, ascii117, ascii118, ascii119, ascii120, ascii121,
	ascii122, ascii123, ascii124, ascii125, ascii126 ]

TEXTALIGNMENT = 'centre' #default value
TEXTANGLE = PI/4 #default value
TEXTWIDTH = 1.0 #default value
TEXTHEIGHT = 1.0 #default value
TEXTSPACING = 0.25 #default value
FONTWIDTH = 4.0 #default value
FONTHEIGHT = 8.0 #default value
FONTSPACING = 1.0 #default value

def charpols (charlist):
    return [MyFont[k] for k in map(lambda x: ord(x)-32, charlist)]

TEXT = COMP([ STRUCT, CAT, DISTR, CONS([ charpols, K(T(1)(FONTSPACING + FONTWIDTH)) ]), CHARSEQ ])

def TEXTWITHATTRIBUTES (TEXTALIGNMENT='centre', TEXTANGLE=0, TEXTWIDTH=1.0, TEXTHEIGHT=2.0, TEXTSPACING=0.25):
    PRINT(TEXTALIGNMENT)
    ALIGN = IF([ K(TEXTALIGNMENT == 'centre'),
				 COMP([ APPLY, CONS([ COMP([ T(1), RAISE(DIV), CONS([ SIZE(1), K(-2) ]) ]), ID ]) ]),
				 IF([ K(TEXTALIGNMENT == 'right'), 
					  COMP([ APPLY, CONS([ COMP([ T(1), RAISE(DIFF), SIZE(1) ]), ID ]) ]), 
					  ID ]) ])
    HANDLE = CONS([
				COMP([ AA(S([1,2])([TEXTWIDTH/FONTWIDTH, TEXTHEIGHT/FONTHEIGHT])), charpols ]),
				K(T(1)(TEXTSPACING + TEXTWIDTH)) ])
    fn = COMP([ R([1,2])(TEXTANGLE), ALIGN, STRUCT, CAT, DISTR, HANDLE, CHARSEQ ])
    return fn

if __name__=="__main__":
    VIEW(TEXT('PLASM'))
    VIEW(TEXTWITHATTRIBUTES('centre',PI/4,0.1,0.4,0.05)('PLASM'))
