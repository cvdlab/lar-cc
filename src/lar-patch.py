from pyplasm import *
from scipy import *
import os,sys
""" import modules from larcc/lib """
sys.path.insert(0, 'lib/py/')
from input_test import *
from larcc import *


hierarchy = [
 [[4, 5, 28, 41, 63, 77, 86, 99, 102, 105],
  [13, 23, 49, 53, 80, 81, 83, 91, 98],
  [12],
  [77, 48],
  [19],
  [95, 47, 31]],
 [[8, 14, 27, 36, 37, 44, 62, 69, 75, 85, 94],
  [1, 24, 40, 50, 60, 61, 65, 87, 92],
  [25],
  [3],
  [66],
  [10, 11, 52, 89, 108]],
 [[20, 21, 30, 67, 70, 88, 90, 96, 107, 110],
  [0, 51, 58, 59, 68, 78, 84, 93],
  [29],
  [56],
  [34],
  [16, 76]],
 [[7, 9, 17, 18, 35, 38, 39, 46, 64, 79, 100, 101, 103, 111],
  [22, 32, 33, 42, 45, 55, 73, 74, 97, 106, 109, 112],
  [6],
  [43, 57],
  [82],
  [104, 54]],
 [[], [], [2], [], [15, 26]]
 ]

def to_pair (item):
  if type(item) is int:
    origin = tuple(AA(min)(TRANS([V[v] for v in FV[item]])))
    pair = (origin, [item])
    return pair
  else:
    items = [x for x in item if type(x) is not list or len(x) > 0]
    pairs = AA(to_pair)(items)
    origins = AA(S1)(pairs)
    origin = (min(AA(S1)(origins)), min(AA(S2)(origins)))
    pair = (origin, pairs)
    return pair
    

pair = to_pair(hierarchy)

print pair

def to_struct (params):
  p_origin, parent = params
  c_origin, c_pairs = parent
  origin = VECTDIFF([c_origin, p_origin])
  if type(c_pairs[0]) is int:
    model = ([V[v] for v in FV[c_pairs[0]]], [range(len(FV[c_pairs[0]]))]) ## ADD EDGES HERE
    return Struct([t(*origin), model])
  else:
    return Struct(  [t(*c_origin)] + AA(to_struct)(DISTL([c_origin, c_pairs]))    )


output = to_struct([pair[0], pair])

VIEW(SKEL_1(STRUCT(MKPOLS( struct2lar(output) ))))


