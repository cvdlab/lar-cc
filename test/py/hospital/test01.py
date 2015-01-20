import sys
sys.path.insert(0, 'lib/py/')
from hospital import *

V,FV,EV = struct2lar(groundFloor)
VIEW(STRUCT(MKPOLS((V,EV))))

V,FV,EV = struct2lar(groundFloor,metric)
VIEW(STRUCT(MKPOLS((V,EV))))
