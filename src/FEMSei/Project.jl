function Project(backend,FTB,Fe::ScalarElement,Grid,QuadOrd,Jacobi,F)
  QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
  Weights = QQ.Weights
  Points = QQ.Points
  fRef  = zeros(Fe.Comp,Fe.DoF,length(Weights))

  p=zeros(Fe.NumG)
  for i = 1 : length(Weights)
    for iComp = 1 : Fe.Comp
      for iD = 1 : Fe.DoF
        fRef[iComp,iD,i] = Fe.phi[iD,iComp](Points[i,1],Points[i,2])
      end
    end
  end
  for iF = 1 : Grid.NumFaces
    pLoc = zeros(Fe.DoF)
    for i = 1 : length(Weights)
      _, detJ, X = Jacobi(Grid.Type,QQ.Points[i,1],QQ.Points[i,2],Grid.Faces[iF], Grid)
      pLoc += abs(detJ)*Weights[i]*(fRef[:,:,i]*F(X[1],X[2],X[3]))
    end
    @. p[Fe.Glob[:,iF]] += pLoc
  end
  p = Fe.M\p
  return p
end
