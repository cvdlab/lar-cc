from larlib import *

filename = "test/svg/boolean/shells.lines"
lines = lines2lines(filename)

V,FV,EV,polygons = larFromLines(lines,False)
VIEW(STRUCT(MKPOLS((V,EV))))