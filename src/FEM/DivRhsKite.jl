function Div1Rhs!(backend,FTB,Div,FeT::ScalarElement,u,uFeF::HDivKiteDElement,
  Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  DivuFRef  = zeros(FeT.Comp,uFeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iDoF = 1 : FeT.DoF
        fTRef[iComp,iDoF,iQ] = FeT.phi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iDoF = 1 : uFeF.DoF
        DivuFRef[iComp,iDoF,iQ] = uFeF.Divphi[iDoF,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  fTRefXP  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  fTRefYP  = zeros(FeT.Comp,FeT.DoF,NumQuadL)
  uFRefXP  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  uFRefYP  = zeros(FeT.Comp,uFeF.DoF,NumQuadL)
  @inbounds for iQ = 1 : NumQuadL
#   @inbounds for iComp = 1 : uFeF.Comp
      @inbounds for iDoF = 1 : uFeF.DoF
        uFRefXP[1,iDoF,iQ] = uFeF.phi[iDoF,1](1.0,PointsL[iQ])
        uFRefYP[1,iDoF,iQ] = uFeF.phi[iDoF,2](PointsL[iQ],1.0)
      end
#   end
    @inbounds for iDoF = 1 : FeT.DoF
      fTRefXP[1,iDoF,iQ] = FeT.phi[iDoF,1](1.0,PointsL[iQ])
      fTRefYP[1,iDoF,iQ] = FeT.phi[iDoF,1](PointsL[iQ],1.0)
    end
  end
  DivLoc = zeros(FeT.DoF)
  uF = zeros(uFeF.DoF)

  @inbounds for iF = 1 : Grid.NumFaces
    DivLoc .= 0
    @inbounds for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uF[iDoF] = u[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      DivuLoc = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF 
        DivuLoc += DivuFRef[1,iDoF,iQ] * uF[iDoF]
      end  
      @inbounds for iDoF = 1 : FeT.DoF
        DivLoc[iDoF] -= Weights[iQ] * fTRef[1,iDoF,iQ] * DivuLoc 
      end
    end
    @inbounds for iQ = 1 : NumQuadL
      uLocX = 0.0 
      uLocY = 0.0 
      @inbounds for iDoF = 1 : uFeF.DoF
        uLocX += uFRefXP[1,iDoF,iQ] * uF[iDoF]
        uLocY += uFRefYP[1,iDoF,iQ] * uF[iDoF]
      end
      @inbounds for iDoF = 1 : FeT.DoF
        DivLoc[iDoF] += WeightsL[iQ] * (fTRefXP[1,iDoF,iQ] * uLocX +
        fTRefYP[1,iDoF,iQ] * uLocY)
      end
    end   
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]  
      Div[ind] += DivLoc[iDoF]
    end  
  end
end
