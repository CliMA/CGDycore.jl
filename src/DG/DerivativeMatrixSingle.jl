function DerivativeMatrixSingle(OrdPoly)
  if OrdPoly == 0
    xw = zeros(1)
    wmat = 2 * ones(1)
  else    
    (xw,wmat)=gausslobatto(OrdPoly+1)
  end
  DS = zeros(OrdPoly+1,OrdPoly+1)
  for i = 1 : OrdPoly + 1
    for j = 1 : OrdPoly + 1
      DS[i,j] = DLagrange(xw[i],xw,j)
      if abs(DS[i,j]) <= 1.e-12
        DS[i,j] = 0.0
      end  
    end
  end
  w = vec(wmat)
  DW = -inv(diagm(w)) * DS' * diagm(w)
  DV = 2 * DS
  DV[1,1] = 2 * DS[1,1] + 1 / w[1]
  DV[OrdPoly+1,OrdPoly+1] = 2 * DS[OrdPoly+1,OrdPoly+1] - 1 / w[OrdPoly+1]
  DVT=DV'
  B = zeros(OrdPoly+1,OrdPoly+1)
  B[1,1] = -1 / w[1]
  B[OrdPoly+1,OrdPoly+1] = 1 / w[OrdPoly+1]
  return (DW,DS,DV,DVT,B)
end
