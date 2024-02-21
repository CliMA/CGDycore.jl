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
    pLoc = zeros(Fe.Comp,Fe.DoF)
    for i = 1 : length(Weights)
      _, detJ, X = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      pLoc += abs(detJ)*Weights[i]*(fRef[:,:,i]*F(X[1],X[2],X[3]))
    end
    @. p[Fe.Glob[:,iF]] += pLoc[Fe.Comp,:]
  end
  p = Fe.M\p
  return p
end

function Project(backend,FTB,Fe::HDivKiteDElement,Grid,QuadOrd,Jacobi,F)
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
      DF, detJ, X = Jacobi(Grid.Type,Points[i,1],Points[i,2],Grid.Faces[iF], Grid)
      pLoc += sign(detJ)*Weights[i]*(fRef[:,:,i]' * (DF' * F(X[1],X[2],X[3])))
    end
    @. p[Fe.Glob[:,iF]] += pLoc[:]
  end
  p1 = Fe.M\p
  return p,p1
end


function ComputeScalar(backend,FTB,Fe::ScalarElement,Grid,p)
  fRef  = zeros(Fe.Comp,Fe.DoF)

  for iComp = 1 : Fe.Comp
    for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](0.0,0.0)
    end
  end
  pM = zeros(Grid.NumFaces,1)
  for iF = 1 : Grid.NumFaces
    pLoc = p[Fe.Glob[:,iF]]
    pM[iF,:] = fRef[:,:]*pLoc
  end
  return pM
end

function ComputeVector(backend,FTB,Fe::HDivKiteDElement,Grid,Jacobi,p)
  fRef  = zeros(Fe.Comp,Fe.DoF)

  for iComp = 1 : Fe.Comp
    for iD = 1 : Fe.DoF
      fRef[iComp,iD] = Fe.phi[iD,iComp](0.0,0.0)
    end
  end

  pM = zeros(Grid.NumFaces,3)
  for iF = 1 : Grid.NumFaces
    DF, detJ, X = Jacobi(Grid.Type,0.0,0.0,Grid.Faces[iF], Grid)
    pLoc = p[Fe.Glob[:,iF]]
    pM[iF,:] = 1/detJ * DF * (fRef[:,:]*pLoc)
  end
  return pM
end
