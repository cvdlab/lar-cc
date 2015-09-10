"""Automatic construction of a simplified 3D building from 2D layout"""
import sys
PATH = "/Users/paoluzzi/Documents/RICERCA/sogei/edifici/"
sys.path.insert(0, PATH)

from buildings import *
from iot3d import *

# LAR models (absolute coordinates)
ala_est = larEmbed(1)(polyline2lar(rects2polylines(eastRooms) + eastTip))
ala_sud = larEmbed(1)(polyline2lar(rects2polylines(southRooms) + southTip))
ala_ovest = larEmbed(1)(polyline2lar(rects2polylines(westRooms) + westTip))
ala_nord = larEmbed(1)(polyline2lar(rects2polylines(northRooms) + northTip)) 
ascensori = larEmbed(1)(polyline2lar(elevators))
spazioComune = larEmbed(1)(polyline2lar(AA(REVERSE)(newLanding)))

# test of input consistency (flat assembly of LAR models)
pianoTipo = Struct([ala_est,ala_sud,ala_ovest,ala_nord,ascensori,spazioComune],"pianoTipo")
##VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(struct2lar(pianoTipo))))
##VIEW(SKEL_1(EXPLODE(1.2,1.2,1.2)(MKPOLS(struct2lar(pianoTipo)))))

# LAR to structs
Ala_est = Struct(lar2Structs(ala_est),"Ala_est")
Ala_sud = Struct(lar2Structs(ala_sud),"Ala_sud")
Ala_ovest = Struct(lar2Structs(ala_ovest),"Ala_ovest")
Ala_nord = Struct(lar2Structs(ala_nord),"Ala_nord")
Ascensori = Struct(lar2Structs(ascensori),"Ascensori")
SpazioComune = Struct(lar2Structs(spazioComune),"SpazioComune")

model = struct2lar(Ala_est)
##VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(model)))
for k,room in enumerate(Ala_est.body):
    print room.body

# hierarchical assembly of simplest LAR models
TreAli = Struct([Ala_est,Ala_sud,Ala_ovest,Ascensori,SpazioComune],"TreAli")
##VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(struct2lar(TreAli))))
Ali_B_C_D = Struct(6*[TreAli, t(0,0,30)],"Ali_B_C_D")
##VIEW(EXPLODE(1.2,1.2,1.2)(MKPOLS(struct2lar(Ali_B_C_D))))

Ala_A = Struct(4*[Ala_nord,  t(0,0,30)],"Ala_A")
ALA_A = EXPLODE(1.2,1.2,1.2)(MKPOLS(struct2lar(Ala_A)))
##VIEW(EXPLODE(1.2,1.2,1.2)([ ALA_A, COLOR(BLUE)(SKEL_1(ALA_A)) ]))

Edificio = Struct([ Ala_A, Ali_B_C_D ],"Edificio")
out = struct2lar(Edificio)
if len(out)==2: V,FV = out
else: V,FV,EV = out

EDIFICIO = STRUCT(MKPOLS((V,FV)))
##VIEW(STRUCT([ EDIFICIO, COLOR(BLUE)(SKEL_1(EDIFICIO)) ]))

VV = AA(LIST)(range(len(V)))
SK = SKEL_1(EDIFICIO)
##VIEW(larModelNumbering(1,1,1)(V,[VV,[],FV],STRUCT([ COLOR(BLUE)(SK), EDIFICIO ]),20))

Va,EVa = struct2lar(TreAli)
Vb,EVb = struct2lar(Ala_A)
a = PROD([ SKEL_1(STRUCT(MKPOLS( ([ v[:2] for v in Va], EVa) ))), QUOTE(5*[30]) ])
b = PROD([ SKEL_1(STRUCT(MKPOLS( ([ v[:2] for v in Vb], EVb) ))), QUOTE(3*[30]) ])
glass = MATERIAL([1,0,0,0.3,  0,1,0,0.3,  0,0,1,0.3, 0,0,0,0.3, 100])
##VIEW(glass(STRUCT([a,b])))

##VIEW(STRUCT([ glass(STRUCT([a,b])), EDIFICIO, COLOR(BLUE)(SKEL_1(EDIFICIO)) ]))


printStruct2GeoJson(PATH,pianoTipo)
printStruct2GeoJson(PATH,Edificio)
