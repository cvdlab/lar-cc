def plan2mockup(model,h,n):
   V,FV = model
   EV = face2edge(FV)
   eE,iP = bUnit_to_eEiP(FV,EV)
   modEe1D = V, [EV[e] for e in eE]
   modIp1D = V, [EV[e] for e in iP]
   modVert0D = larExtrude1( VOID, 6*[1] )
   modVert1D = larExtrude1( VOID, 6*[1] )
   horParts = larModelProduct([ model, modVert0D ])
   
VIEW(EXPLODE(1.2,1.2,1)(MKPOLS((X,FX))))
