from larlib import *
from meshpy.tet import MeshInfo, build
    
def faces(tet):
	v1,v2,v3,v4 = tet
	return AA(sorted)([(v1,v2,v3),(v2,v3,v4),(v3,v4,v1),(v4,v1,v2)])

def edges(tria):
	v1,v2,v3 = tria
	return AA(sorted)([(v1,v2),(v1,v3),(v2,v3)])

def brep2lar(larBrep):
    V,FV = larBrep    
    mesh_info = MeshInfo()
    mesh_info.set_points(V)
    mesh_info.set_facets(FV)
    mesh = build(mesh_info)
    W = [v for h,v in enumerate(mesh.points)]
    CW = [tet for k,tet in enumerate(mesh.elements)]
    FW = sorted(set(AA(tuple)(CAT(AA(faces)(CW)))))
    EW = sorted(set(AA(tuple)(CAT(AA(edges)(FW)))))
    return W,CW,FW,EW

V = [(0,0,0), (2,0,0), (2,2,0), (0,2,0), (0,0,12), (2,0,12), (2,2,12), (0,2,12)]
FV = [[0,1,2,3],[4,5,6,7], [0,4,5,1], [1,5,6,2], [2,6,7,3], [3,7,4,0]]

VV = AA(LIST)(range(len(V)))
hpc = SKEL_1(STRUCT(MKPOLS((V,FV))))
VIEW(larModelNumbering(1,1,1)(V,[VV,FV],hpc,4))

W,CW,FW,EW = brep2lar((V,FV))

WW = AA(LIST)(range(len(W)))
hpc = SKEL_1(STRUCT(MKPOLS((W,FW))))
VIEW(larModelNumbering(1,1,1)(W,[WW,FW,EW],hpc,2))

    
V,(VV,EV,FV,CV) = larCuboids((1,1,2),True)
VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS((V,FV))))