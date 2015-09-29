""" Test of pruning clusters of close vertices """
from larlib import *
from scipy import rand
from scipy.spatial import cKDTree
POINTS = 1000
RADIUS = 0.01

pts = [rand(2).tolist() for k in range(POINTS)]
VIEW(STRUCT(AA(MK)(pts)))
V,close,clusters,vmap = pruneVertices(pts,RADIUS)
circles = [T([1,2])(pts[h])(CIRCUMFERENCE(RADIUS)(18)) for h,k in close]
convexes = [JOIN(AA(MK)([pts[v] for v in cluster])) for cluster in clusters]
W = COLOR(CYAN)(STRUCT(AA(MK)(V)))
VIEW(STRUCT(AA(MK)(pts)+AA(COLOR(YELLOW))(circles)))
VIEW(STRUCT(AA(COLOR(RED))(convexes)+AA(MK)(pts)+AA(COLOR(YELLOW))(circles)+[W]))
