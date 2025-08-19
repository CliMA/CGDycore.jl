function GradKinHeightKite!(backend,FTB,Rhs,h,hFeF::ScalarElement,u,uFeF::HDivKiteDElement,
  FeT::HDivKiteDElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  NumQuad, Weights, Points = QuadRule(Grid.Type,QuadOrd)
  fFRef  = zeros(hFeF.Comp,hFeF.DoF,length(Weights))
  fTRef  = zeros(hFeF.Comp,FeT.DoF,length(Weights))
  uFRef  = zeros(uFeF.Comp,uFeF.DoF,length(Weights))

  for iQ = 1 : length(Weights)
    for iComp = 1 : hFeF.Comp
      for iDoF = 1 : FeT.DoF
        fTRef[iComp,iDoF,iQ] = FeT.Divphi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : hFeF.Comp
      for iDoF = 1 : hFeF.DoF
        fFRef[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : uFeF.Comp
      for iDoF = 1 : uFeF.DoF
        uFRef[iComp,iDoF,iQ] = uFeF.phi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fTRefX  = zeros(hFeF.Comp,FeT.DoF,length(WeightsL))
  fTRefY  = zeros(hFeF.Comp,FeT.DoF,length(WeightsL))
  hFRefX  = zeros(hFeF.Comp,hFeF.DoF,length(WeightsL))
  hFRefY  = zeros(hFeF.Comp,hFeF.DoF,length(WeightsL))
  uFRefX  = zeros(uFeF.Comp,uFeF.DoF,length(WeightsL))
  uFRefY  = zeros(uFeF.Comp,uFeF.DoF,length(WeightsL))
  for iQ = 1 : length(WeightsL)
    for iComp = 1 : hFeF.Comp
      for iDoF = 1 : hFeF.DoF
        hFRefX[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](-1.0,PointsL[iQ])
        hFRefY[iComp,iDoF,iQ] = hFeF.phi[iDoF,iComp](PointsL[iQ],-1.0)
      end
    end
    for iDoF = 1 : uFeF.DoF
      uFRefX[1,iDoF,iQ] = uFeF.phi[iDoF,1](-1.0,PointsL[iQ])
      uFRefX[2,iDoF,iQ] = uFeF.phi[iDoF,2](-1.0,PointsL[iQ])
      uFRefY[1,iDoF,iQ] = uFeF.phi[iDoF,1](PointsL[iQ],-1.0)
      uFRefY[2,iDoF,iQ] = uFeF.phi[iDoF,2](PointsL[iQ],-1.0)
    end
    for iDoF = 1 : FeT.DoF
      fTRefX[1,iDoF,iQ] = FeT.phi[iDoF,1](-1.0,PointsL[iQ])
      fTRefY[1,iDoF,iQ] = FeT.phi[iDoF,2](PointsL[iQ],-1.0)
    end
  end
  GradLoc = zeros(FeT.DoF)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  for iF = 1 : Grid.NumFaces
    GradLoc .= 0
    for iQ = 1 : length(Weights)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1 = 0.0
      u2 = 0.0
      for iDoF = 1 : uFeF.DoF
        ind = uFeF.Glob[iDoF,iF]
        @views u1 += uFRef[1,iDoF,iQ] * u[ind]
        @views u2 += uFRef[2,iDoF,iQ] * u[ind]
      end
      uLoc1 = 1 / detDFLoc * (DF[1,1] * u1 + DF[1,2] * u2)
      uLoc2 = 1 / detDFLoc * (DF[2,1] * u1 + DF[2,2] * u2)
      uLoc3 = 1 / detDFLoc * (DF[3,1] * u1 + DF[3,2] * u2)
      hLoc = 0.5 * (uLoc1 * uLoc1 + uLoc2 * uLoc2 + uLoc3 * uLoc3)
      for iDoF = 1 : hFeF.DoF 
        ind = hFeF.Glob[iDoF,iF]
        hLoc += fFRef[1,iDoF,iQ] * h[ind]
      end  
      GradLoc += sign(detDFLoc) * Weights[iQ] * fTRef[1,:,iQ] * hLoc
    end
    for iQ = 1 : length(WeightsL)
      hLocX = 0.0 
      hLocY = 0.0 
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,-1.0,PointsL[iQ],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1X = 0.0
      u2X = 0.0
      for iDoF = 1 : uFeF.DoF
        ind = uFeF.Glob[iDoF,iF]
        u1X += uFRefX[1,iDoF,iQ] * u[ind]
        u2X += uFRefX[2,iDoF,iQ] * u[ind]
      end
      uLoc1X = 1 / detDFLoc * (DF[1,1] * u1X + DF[1,2] * u2X)
      uLoc2X = 1 / detDFLoc * (DF[2,1] * u1X + DF[2,2] * u2X)
      uLoc3X = 1 / detDFLoc * (DF[3,1] * u1X + DF[3,2] * u2X)
      hLocX = 0.5 * (uLoc1X * uLoc1X + uLoc2X * uLoc2X + uLoc3X * uLoc3X)
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,PointsL[iQ],-1.0,Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      u1Y = 0.0
      u2Y = 0.0
      for iDoF = 1 : uFeF.DoF
        ind = uFeF.Glob[iDoF,iF]
        u1Y += uFRefY[1,iDoF,iQ] * u[ind]
        u2Y += uFRefY[2,iDoF,iQ] * u[ind]
      end
      uLoc1Y = 1 / detDFLoc * (DF[1,1] * u1Y + DF[1,2] * u2Y)
      uLoc2Y = 1 / detDFLoc * (DF[2,1] * u1Y + DF[2,2] * u2Y)
      uLoc3Y = 1 / detDFLoc * (DF[3,1] * u1Y + DF[3,2] * u2Y)
      hLocY = 0.5 * (uLoc1Y * uLoc1Y + uLoc2Y * uLoc2Y + uLoc3Y * uLoc3Y)
      for iDoF = 1 : hFeF.DoF
        ind = hFeF.Glob[iDoF,iF]
        hLocX += hFRefX[1,iDoF,iQ] * h[ind]
        hLocY += hFRefY[1,iDoF,iQ] * h[ind]
      end
      GradLoc += WeightsL[iQ] * (fTRefX[1,:,iQ] * hLocX + fTRefY[1,:,iQ] * hLocY)
    end   
    for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]  
      Rhs[ind] += GradLoc[iDoF]
    end  
  end
end
