""" make the model of a layout floor """
from larlib import *

""" make the model of a layout floor """

scaleFactor = 83.333

filename = "test/svg/inters/plan.svg"
larModel = svg2lar(filename)
larModel = larApply(s(scaleFactor,scaleFactor))(larModel)
V,FV,EV = larModel
FV[2] += FV[71]      # for now :o)

""" Assembling floor layout generation """
""" Visualization of cell numbering in a 2D complex """
VV = AA(LIST)(range(len(V)))
submodel = STRUCT(MKPOLS((V,EV)))
#VIEW(larModelNumbering(1,1,1)(V,[VV,EV,FV[:-1]],submodel,2.5))

VV = AA(LIST)(range(len(V)))
FE = crossRelation(FV,EV,VV)
chainsToStruct = chain2structs(V,FV,EV,FE)

""" Ala nord """
boxes = [0 for k in range(64)]
point = [0 for k in range(64)]
boxes[0] = array([[0.431, 0.607], [0.474, 0.91]])*scaleFactor #[V[k] for k in [39,208]]
boxes[1] = array([[0.416, 0.657], [0.372, 0.953]])*scaleFactor #[V[k] for k in [162,39]]
boxes[2] = array([[0.416, 0.627], [0.431, 0.986]])*scaleFactor #[V[k] for k in [206,247]]
boxes[3] = array([[0.431, 0.607], [0.448, 0.627]])*scaleFactor #[V[k] for k in [39,7]]
boxes[4] = array([[0.431, 0.91], [0.494, 0.929]])*scaleFactor  #[V[k] for k in [213,234]]
boxes[5] = array([[0.431, 0.97], [0.466, 1.0]])*scaleFactor #[V[k] for k in [58,88]]
boxes[27] = array([[0.416, 0.627], [0.372, 0.657]])*scaleFactor #[V[k] for k in [110,82]]

point[0] = array([0.394, 0.9625])*scaleFactor #CCOMB([V[k] for k in [190,197]])
point[1] = array([0.4525, 0.9325])*scaleFactor #CCOMB([V[k] for k in [166,159]])

piano1_superficieUtile_zonaNord_uffici_destra = subComplexInBox(V,FV,EV,boxes[0])[1]
piano1_superficieUtile_zonaNord_uffici_sinistra = subComplexInBox(V,FV,EV,boxes[1])[1]
piano1_connettivo_orizzontale_zonaNord = subComplexInBox(V,FV,EV,boxes[2])[1]
piano1_connettivo_verticale_zonaNord_ascensore = subComplexInBox(V,FV,EV,boxes[3])[1]
piano1_connettivo_verticale_zonaNord_ascensore += subComplexInBox(V,FV,EV,boxes[4])[1]
piano1_connettivo_verticale_zonaNord_scale = subComplexInBox(V,FV,EV,boxes[5])[1]
piano1_superficieUtile_zonaNord_servizi = subComplexAroundPoint(V,FV,EV,FE,point[0])[1]
piano1_superficieUtile_zonaNord_servizi += subComplexAroundPoint(V,FV,EV,FE,point[1])[1]
piano1_superficieUtile_zonaNord_servizi += subComplexInBox(V,FV,EV,boxes[27])[1]

piano1N = [piano1_superficieUtile_zonaNord_uffici_destra, piano1_superficieUtile_zonaNord_uffici_sinistra, piano1_connettivo_orizzontale_zonaNord, piano1_connettivo_verticale_zonaNord_ascensore, piano1_connettivo_verticale_zonaNord_scale, piano1_superficieUtile_zonaNord_servizi]
""" Ala nord """
piano1N_nomi = ["piano1_superficieUtile_zonaNord_uffici_destra", "piano1_superficieUtile_zonaNord_uffici_sinistra", "piano1_connettivo_orizzontale_zonaNord", "piano1_connettivo_verticale_zonaNord_ascensore", "piano1_connettivo_verticale_zonaNord_scale", "piano1_superficieUtile_zonaNord_servizi"]
piano1N_categorie = ["uffici","uffici","corridoi","ascensori","scale","servizi"]
p1N = zip(piano1N,piano1N_nomi,piano1N_categorie)
piano1_zonaNord = Struct(AA(chainsToStruct)(p1N),"piano1_zonaNord","ala")
#VIEW(SKEL_1(STRUCT(MKPOLS(struct2lar(piano1_zonaNord)))))
    
nord = CAT([cells2hpcs(V,FV,chain,k) for k,chain in enumerate(piano1N)])
#VIEW(EXPLODE(1.2,1.2,1.2)(nord))

""" Ala est """
boxes[6] = array([[0.019, 0.533], [0.376, 0.577]])*scaleFactor #[V[k] for k in [241,29]]
boxes[7] = array([[0.07, 0.474], [0.343, 0.518]])*scaleFactor #[V[k] for k in [264,148]]
boxes[8] = array([[0.013, 0.518], [0.376, 0.533]])*scaleFactor #[V[k] for k in [22,63]]
boxes[9] = array([[0.376, 0.533], [0.39, 0.549]])*scaleFactor #[V[k] for k in [63,92]]
boxes[10] = array([[0.001, 0.474], [0.07, 0.518]])*scaleFactor #[V[k] for k in [263,265]]
boxes[11] = array([[0.343, 0.474], [0.376, 0.518]])*scaleFactor #[V[k] for k in [84,149]]

point[2] = array([0.015, 0.5535])*scaleFactor #CCOMB([V[k] for k in [228,14]])

piano1_superficieUtile_zonaEst_uffici_destra = subComplexInBox(V,FV,EV,boxes[6])[1]
piano1_superficieUtile_zonaEst_uffici_sinistra = subComplexInBox(V,FV,EV,boxes[7])[1]
piano1_connettivo_orizzontale_zonaEst = subComplexInBox(V,FV,EV,boxes[8])[1]
piano1_connettivo_verticale_zonaEst_ascensore = subComplexInBox(V,FV,EV,boxes[9])[1]
piano1_connettivo_verticale_zonaEst_scale = subComplexAroundPoint(V,FV,EV,FE,point[2])[1]
piano1_superficieUtile_zonaEst_servizi = subComplexInBox(V,FV,EV,boxes[10])[1]
piano1_superficieUtile_zonaEst_servizi += subComplexInBox(V,FV,EV,boxes[11])[1]

piano1E = [piano1_superficieUtile_zonaEst_uffici_destra, piano1_superficieUtile_zonaEst_uffici_sinistra, piano1_connettivo_orizzontale_zonaEst, piano1_connettivo_verticale_zonaEst_ascensore, piano1_connettivo_verticale_zonaEst_scale, piano1_superficieUtile_zonaEst_servizi]
""" Ala est """
piano1E_nomi = ["piano1_superficieUtile_zonaEst_uffici_destra", "piano1_superficieUtile_zonaEst_uffici_sinistra", "piano1_connettivo_orizzontale_zonaEst", "piano1_connettivo_verticale_zonaEst_ascensore", "piano1_connettivo_verticale_zonaEst_scale", "piano1_superficieUtile_zonaEst_servizi"]
piano1E_categorie = ["uffici","uffici","corridoi","ascensori","scale","servizi"]
p1E = zip(piano1E,piano1E_nomi, piano1E_categorie)
piano1_zonaEst = Struct(AA(chainsToStruct)(p1E), "piano1_zonaEst", "ala")
#VIEW(SKEL_1(STRUCT(MKPOLS(struct2lar(piano1_zonaEst)))))

est = CAT([cells2hpcs(V,FV,chain,k) for k,chain in enumerate(piano1E)])
#VIEW(EXPLODE(1.2,1.2,1.2)(est + nord))

""" Ala sud """
boxes[12] = array([[0.467, 0.138], [0.423, 0.476]])*scaleFactor #[V[k] for k in [252,47]]
boxes[13] = array([[0.482, 0.145], [0.525, 0.445]])*scaleFactor #[V[k] for k in [241,126]]
boxes[14] = array([[0.482, 0.476], [0.467, 0.116]])*scaleFactor #[V[k] for k in [254,232]]
boxes[15] = array([[0.449, 0.476], [0.467, 0.493]])*scaleFactor #[V[k] for k in [40,237]]
boxes[16] = array([[0.431, 0.101], [0.467, 0.131]])*scaleFactor #[V[k] for k in [259,2]]
boxes[17] = array([[0.482, 0.445], [0.525, 0.476]])*scaleFactor #[V[k] for k in [155,248]]
boxes[18] = array([[0.525, 0.104], [0.482, 0.145]])*scaleFactor #[V[k] for k in [111,241]]

piano1_superficieUtile_zonaSud_uffici_destra = subComplexInBox(V,FV,EV,boxes[12])[1]
piano1_superficieUtile_zonaSud_uffici_sinistra = subComplexInBox(V,FV,EV,boxes[13])[1]
piano1_connettivo_orizzontale_zonaSud = subComplexInBox(V,FV,EV,boxes[14])[1]
piano1_connettivo_verticale_zonaSud_ascensore = subComplexInBox(V,FV,EV,boxes[15])[1]
piano1_connettivo_verticale_zonaSud_scale = subComplexInBox(V,FV,EV,boxes[16])[1]
piano1_superficieUtile_zonaSud_servizi = subComplexInBox(V,FV,EV,boxes[17])[1]
piano1_superficieUtile_zonaSud_servizi += subComplexInBox(V,FV,EV,boxes[18])[1]

piano1S = [piano1_superficieUtile_zonaSud_uffici_destra, piano1_superficieUtile_zonaSud_uffici_sinistra, piano1_connettivo_orizzontale_zonaSud, piano1_connettivo_verticale_zonaSud_ascensore, piano1_connettivo_verticale_zonaSud_scale, piano1_superficieUtile_zonaSud_servizi]
""" Ala sud """
piano1S_nomi = ["piano1_superficieUtile_zonaSud_uffici_destra", "piano1_superficieUtile_zonaSud_uffici_sinistra", "piano1_connettivo_orizzontale_zonaSud", "piano1_connettivo_verticale_zonaSud_ascensore", "piano1_connettivo_verticale_zonaSud_scale", "piano1_superficieUtile_zonaSud_servizi"]
piano1S_categorie = ["uffici","uffici","corridoi","ascensori","scale","servizi"]
p1S = zip(piano1S,piano1S_nomi, piano1S_categorie)
piano1_zonaSud = Struct(AA(chainsToStruct)(p1S), "piano1_zonaSud", "ala")
#VIEW(SKEL_1(STRUCT(MKPOLS(struct2lar(piano1_zonaSud)))))
    
sud = CAT([cells2hpcs(V,FV,chain,k) for k,chain in enumerate(piano1S)])
#VIEW(EXPLODE(1.2,1.2,1.2)(est + nord + sud))

""" Ala ovest """
boxes[19] = array([[0.521, 0.526], [0.963, 0.568]])*scaleFactor #[V[k] for k in [169,202]]
boxes[20] = array([[0.555, 0.584], [0.955, 0.627]])*scaleFactor #[V[k] for k in [12,23]]
boxes[21] = array([[0.521, 0.568], [0.985, 0.584]])*scaleFactor #[V[k] for k in [209,204]]
boxes[22] = array([[0.506, 0.551], [0.521, 0.568]])*scaleFactor #[V[k] for k in [89,209]]
boxes[23] = array([[0.808, 0.504], [0.828, 0.526]])*scaleFactor #[V[k] for k in [270,77]]
boxes[24] = array([[0.955, 0.584], [0.997, 0.627]])*scaleFactor #[V[k] for k in [220,24]]
boxes[25] = array([[0.521, 0.584], [0.555, 0.627]])*scaleFactor #[V[k] for k in [11,144]]
boxes[26] = array([[1.0, 0.533], [0.97, 0.568]])*scaleFactor #[V[k] for k in [233,201]]

piano1_superficieUtile_zonaOvest_uffici_destra = subComplexInBox(V,FV,EV,boxes[19])[1]
piano1_superficieUtile_zonaOvest_uffici_sinistra = subComplexInBox(V,FV,EV,boxes[20])[1]
piano1_connettivo_orizzontale_zonaOvest = subComplexInBox(V,FV,EV,boxes[21])[1]
piano1_connettivo_verticale_zonaOvest_ascensore = subComplexInBox(V,FV,EV,boxes[22])[1]
piano1_connettivo_verticale_zonaOvest_ascensore += subComplexInBox(V,FV,EV,boxes[23])[1]
piano1_superficieUtile_zonaOvest_servizi = subComplexInBox(V,FV,EV,boxes[24])[1]
piano1_superficieUtile_zonaOvest_servizi += subComplexInBox(V,FV,EV,boxes[25])[1]
piano1_connettivo_verticale_zonaOvest_scale = subComplexInBox(V,FV,EV,boxes[26])[1]

piano1O = [piano1_superficieUtile_zonaOvest_uffici_destra, piano1_superficieUtile_zonaOvest_uffici_sinistra, piano1_connettivo_orizzontale_zonaOvest, piano1_connettivo_verticale_zonaOvest_ascensore, piano1_connettivo_verticale_zonaOvest_scale, piano1_superficieUtile_zonaOvest_servizi]
""" Ala ovest """
piano1O_nomi = ["piano1_superficieUtile_zonaOvest_uffici_destra", "piano1_superficieUtile_zonaOvest_uffici_sinistra", "piano1_connettivo_orizzontale_zonaOvest", "piano1_connettivo_verticale_zonaOvest_ascensore", "piano1_connettivo_verticale_zonaOvest_scale", "piano1_superficieUtile_zonaOvest_servizi"]
piano1O_categorie = ["uffici","uffici","corridoi","ascensori","scale","servizi"]
p1O = zip(piano1O,piano1O_nomi, piano1O_categorie)
piano1_zonaOvest = Struct(AA(chainsToStruct)(p1O), "piano1_zonaOvest", "ala")
#VIEW(SKEL_1(STRUCT(MKPOLS(struct2lar(piano1_zonaOvest)))))
    
ovest = CAT([cells2hpcs(V,FV,chain,k) for k,chain in enumerate(piano1O)])
#VIEW(EXPLODE(1.2,1.2,1.2)(est + nord + sud + ovest))

""" Centro stella """
piano1_connettivo_orizzontale_centroStella = [2]
piano1_connettivo_verticale_centroStella_scale = [15,26]

piano1C = [[],[],piano1_connettivo_orizzontale_centroStella,[], piano1_connettivo_verticale_centroStella_scale]
centro = CAT([cells2hpcs(V,FV,chain,k) for k,chain in enumerate(piano1C)])
""" Centro stella """
piano1C_nomi = [[],[],"piano1_connettivo_orizzontale_centroStella", [], "piano1_connettivo_verticale_centroStella_scale"]
piano1C_categorie = [[],[],"corridoi",[], "ascensori"]
p1C = zip(piano1C,piano1C_nomi, piano1C_categorie)
piano1_centroStella = Struct(AA(chainsToStruct)(p1C), "piano1_centroStella", "centro")
#VIEW(SKEL_1(STRUCT(MKPOLS(struct2lar(piano1_centroStella)))))

#VIEW(EXPLODE(1.2,1.2,1.2)(est + nord + sud + ovest + centro))
#VIEW(STRUCT(est + nord + sud + ovest + centro))

""" Assemblaggio """
p1 = p1N + p1S + p1E + p1O + p1C

piano1_nomi = ["piano1_zonaNord", "piano1_zonaEst", "piano1_zonaSud", "piano1_zonaOvest", "piano1_centroStella"]
piano1_categorie = ["ala","ala","ala","ala","centro"]
piano1 = Struct(AA(chainsToStruct)(p1), "piano1", "level")

#VIEW(SKEL_1(STRUCT(MKPOLS(struct2lar(piano1)))))

V,boundaryEdges = structBoundaryModel(piano1)
drawing = mkSignedEdges((V,boundaryEdges))
#VIEW(drawing)
        
polylines = boundaryModel2polylines((V,boundaryEdges))
#VIEW(STRUCT(AA(POLYLINE)(polylines)))
    
print boundaryPolylines(piano1)



primoPiano = AA(list)([CAT(AA(S1)(p1N)),CAT(AA(S1)(p1E)),CAT(AA(S1)(p1S)), 
                CAT(AA(S1)(p1O)),CAT(AA(S1)(p1C))])
primoPiano_nomi = ["piano1_zonaNord","piano1_zonaEst","piano1_zonaSud","piano1_zonaOvest","piano1_centroStella"]
primoPiano_categorie = ["ala","ala","ala","ala","centro"]
pianoPrimo = zip(primoPiano, primoPiano_nomi, primoPiano_categorie)
piano_1 = Struct( AA(chainsToStruct)(pianoPrimo), "piano1", "level" )

piano_1_3D = embedStruct(1)(piano_1,"3D")
iot3d.printStruct2GeoJson("./",piano_1_3D)

print "piano_1_3D =",piano_1_3D.category
print "piano_1_3D.body[0] =",piano_1_3D.body[0].category
print "piano_1_3D.body[0].body[0] =",piano_1_3D.body[0].body[0].category
