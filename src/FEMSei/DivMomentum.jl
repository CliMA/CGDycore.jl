#scalar variant Quads
function DivMomentumScalar!(backend,FTB,Rhs,uHDiv,FeHDiv::HDivElement,cDG,FeDG::ScalarElement,FeT::ScalarElement,Grid,
  ElemType::Grids.Quad,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  fuHDiv  = zeros(FeHDiv.DoF,FeHDiv.Comp,NumQuad) 
  fuDG  = zeros(FeDG.DoF,FeDG.DoF,NumQuad) 
  fuT  = zeros(FeT.DoF,FeT.Comp,NumQuad) 
  fuDivHDiv = zeros(FeHDiv.DoF,NumQuad) 
  fuGradDG = zeros(FeDG.DoF,FeDG.DoF,2,NumQuad) 

  #computation of the ansatz functions in the quadrature points
  for iQ = 1 : NumQuad
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fuT[iD,iComp,iQ] = FeT.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeHDiv.DoF
        fuHDiv[iD,iComp,iQ] = FeHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iD = 1 : FeDG.DoF
      fuGradDG[iD,1,1,iQ] = FeDG.Gradphi[iD,1,1](Points[iQ,1],Points[iQ,2])
      fuGradDG[iD,1,2,iQ] = FeDG.Gradphi[iD,1,2](Points[iQ,1],Points[iQ,2])
    end
    for iD = 1 : FeDG.DoF
      fuDG[iD,1,iQ] = FeDG.phi[iD,1](Points[iQ,1],Points[iQ,2])
    end
    for iD = 1 : FeHDiv.DoF
      fuDivHDiv[iD,iQ] = FeHDiv.Divphi[iD,1](Points[iQ,1],Points[iQ,2])
    end
  end

  #local variables
  cDGLoc = zeros(FeDG.DoF)
  uHDivLoc = zeros(FeHDiv.DoF)
  RhsLoc = zeros(FeT.DoF)

  uHDivLocLeft = zeros(FeHDiv.DoF)
  uHDivLocRight = zeros(FeHDiv.DoF)
  cDGLocLeft = zeros(FeDG.DoF)
  cDGLocRight = zeros(FeDG.DoF)

  uuGradDGLoc = zeros(2)
  uuHDivLoc = zeros(2)
 
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
 
  #copy from global variables to local variables  
  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iD = 1 : FeDG.DoF
      ind = FeDG.Glob[iD,iF]  
      cDGLoc[iD] = cDG[ind]
    end 
    @inbounds for iD = 1 : FeHDiv.DoF
      ind = FeHDiv.Glob[iD,iF]  
      uHDivLoc[iD] = uHDiv[ind]
    end  
    @. RhsLoc = 0.0

    for iQ = 1 : NumQuad
      #computation of local variables in a quadrature point
      uuDGLoc = 0.0
      @. uuGradDGLoc = 0.0
      @inbounds for iD = 1 : FeDG.DoF
        uuDGLoc += cDGLoc[iD] * fuDG[iD,1,iQ] 
        @views @. uuGradDGLoc += cDGLoc[iD] * fuGradDG[iD,1,:,iQ] 
      end
      uuDivHDivLoc = 0.0
      @. uuHDivLoc = 0.0
      @inbounds for iD = 1 : FeHDiv.DoF
        uuDivHDivLoc += uHDivLoc[iD] * fuDivHDiv[iD,iQ] 
        uuHDivLoc[1] += uHDivLoc[iD] * fuHDiv[iD,1,iQ]
        uuHDivLoc[2] += uHDivLoc[iD] * fuHDiv[iD,2,iQ] 
      end
      GradDGHDiv = uuGradDGLoc' * uuHDivLoc
      DGDivHDiv = uuDivHDivLoc * uuDGLoc
      innersum = Grid.Faces[iF].Orientation * (DGDivHDiv + GradDGHDiv)
      #product incoming functions and test function
      @inbounds for iD = 1 : FeT.DoF
        product = innersum * fuT[iD,1,iQ]
        RhsLoc[iD] += Weights[iQ] * product
      end 
    end

    #save the new without upwind 
    @inbounds for iD = 1 : FeT.DoF
      ind = FeT.Glob[iD,iF] 
      Rhs[ind] += RhsLoc[iD]
    end  
  end 

  #UPWIND on Edges
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFFRef  = zeros(FeHDiv.DoF, FeHDiv.Comp,NumQuadL,4)
  hFFRef  = zeros(FeDG.DoF,FeDG.DoF,NumQuadL,4)
  fTRef  = zeros(FeT.DoF,FeT.Comp,NumQuadL,4)
  uFFLocL = zeros(FeHDiv.DoF)
  uFFLocR = zeros(FeHDiv.DoF)

  #testfunctions DG scalar 
  for iQ = 1 : NumQuadL 
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iD,iComp,iQ,1] = FeT.phi[iD,iComp](PointsL[iQ],-1.0)
        fTRef[iD,iComp,iQ,2] = FeT.phi[iD,iComp](1.0,PointsL[iQ])
        fTRef[iD,iComp,iQ,3] = FeT.phi[iD,iComp](PointsL[iQ],1.0)
        fTRef[iD,iComp,iQ,4] = FeT.phi[iD,iComp](-1.0,PointsL[iQ])
      end
    end
    #uf RT UHDiv
    for iComp = 1 : FeHDiv.Comp
      for iD = 1 : FeHDiv.DoF
        uFFRef[iD,iComp,iQ,1] = FeHDiv.phi[iD,iComp](PointsL[iQ],-1.0)
        uFFRef[iD,iComp,iQ,2] = FeHDiv.phi[iD,iComp](1.0,PointsL[iQ])
        uFFRef[iD,iComp,iQ,3] = FeHDiv.phi[iD,iComp](PointsL[iQ],1.0)
        uFFRef[iD,iComp,iQ,4] = FeHDiv.phi[iD,iComp](-1.0,PointsL[iQ])
      end
    end
    #hf
    for iComp = 1 : FeDG.Comp
      for iD = 1 : FeDG.DoF
        hFFRef[iD,iComp,iQ,1] = FeDG.phi[iD,iComp](PointsL[iQ],-1.0)
        hFFRef[iD,iComp,iQ,2] = FeDG.phi[iD,iComp](1.0,PointsL[iQ])
        hFFRef[iD,iComp,iQ,3] = FeDG.phi[iD,iComp](PointsL[iQ],1.0)
        hFFRef[iD,iComp,iQ,4] = FeDG.phi[iD,iComp](-1.0,PointsL[iQ])
      end
    end
  end

  #allocation Mloc
  MLoc11 = zeros(FeT.DoF, FeHDiv.DoF)
  MLoc12 = zeros(FeT.DoF, FeHDiv.DoF)
  MLoc21 = zeros(FeT.DoF, FeHDiv.DoF)
  MLoc22 = zeros(FeT.DoF, FeHDiv.DoF)

  @inbounds for iE in 1:Grid.NumEdges
    Edge = Grid.Edges[iE]
    if length(Edge.F) > 1
      iFL = Edge.F[1]
      EdgeTypeL = Edge.FE[1]
      iFR = Edge.F[2]
      EdgeTypeR = Edge.FE[2]
      #computation normales of edges
      @views nBarLocL = Grid.nBar[:, EdgeTypeL]
      @views nBarLocR = Grid.nBar[:, EdgeTypeR]
      #gamma upwind value
      uE = 0.0 
      @inbounds for iD = 1 : FeHDiv.DoF
        ind = FeHDiv.Glob[iD,iFL]  
        uHDivLocLeft[iD] = uHDiv[ind]
        ind = FeHDiv.Glob[iD,iFR]  
        uHDivLocRight[iD] = uHDiv[ind]
      end 
      @inbounds for iD = 1 : FeDG.DoF
        ind = FeDG.Glob[iD,iFL]  
        cDGLocLeft[iD] = cDG[ind]
        ind = FeDG.Glob[iD,iFR]  
        cDGLocRight[iD] = cDG[ind]
      end 
      @inbounds for iQ = 1:length(WeightsL)
        @inbounds for iD = 1 : FeHDiv.DoF
          @views uE += WeightsL[iQ]*(uHDivLocLeft[iD]*(uFFRef[iD,:,iQ,EdgeTypeL]'*nBarLocL))
        end
      end
      #upwind value
      gammaU = 0.5
      gammaLoc = uE > 0 ? gammaU : -gammaU

      @. MLoc11 = 0
      @. MLoc12 = 0
      @. MLoc21 = 0
      @. MLoc22 = 0
     

      @inbounds for iQ in 1:length(WeightsL)
        @inbounds for iD in 1:FeHDiv.DoF
          @views uFFLocL[iD] = uFFRef[iD,:,iQ,EdgeTypeL]' * nBarLocL
          @views uFFLocR[iD] = uFFRef[iD,:,iQ,EdgeTypeR]' * nBarLocR
        end
        hFFL = 0.0
        hFFR = 0.0
        @inbounds for iD in 1:FeDG.DoF
          hFFL += hFFRef[iD,1,iQ,EdgeTypeL] * cDGLocLeft[iD]
          hFFR += hFFRef[iD,1,iQ,EdgeTypeR] * cDGLocRight[iD]
        end
        @inbounds for iDoF  = 1 : FeT.DoF 
          @inbounds for jDoF  = 1 : FeHDiv.DoF
            MLoc11[iDoF,jDoF] += WeightsL[iQ] * (fTRef[iDoF,1,iQ,EdgeTypeL] * hFFL) * uFFLocL[jDoF]'
            MLoc12[iDoF,jDoF] += WeightsL[iQ] * (fTRef[iDoF,1,iQ,EdgeTypeL] * hFFR) * uFFLocR[jDoF]'
            MLoc21[iDoF,jDoF] += WeightsL[iQ] * (fTRef[iDoF,1,iQ,EdgeTypeR] * hFFL) * uFFLocL[jDoF]'
            MLoc22[iDoF,jDoF] += WeightsL[iQ] * (fTRef[iDoF,1,iQ,EdgeTypeR] * hFFR) * uFFLocR[jDoF]'
          end
        end
      end
      #all together over edges
      @inbounds for iD = 1 : FeT.DoF
        indL = FeT.Glob[iD,iFL] 
        @views Rhs[indL] += (+0.5 - gammaLoc) * MLoc11[iD,:]' * uHDivLocLeft
        @views Rhs[indL] += (-0.5 + gammaLoc) * MLoc12[iD,:]' * uHDivLocRight
        indR = FeT.Glob[iD,iFR] 
        @views Rhs[indR] += (+0.5 + gammaLoc) * MLoc21[iD,:]' * uHDivLocLeft
        @views Rhs[indR] += (-0.5 - gammaLoc) * MLoc22[iD,:]' * uHDivLocRight
      end
    end
  end
end

#scalar variant Tri
function DivMomentumScalar!(backend,FTB,Rhs,uHDiv,FeHDiv::HDivElement,cDG,FeDG::ScalarElement,FeT::ScalarElement,Grid,
  ElemType::Grids.Tri,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  fuHDiv  = zeros(FeHDiv.DoF,FeHDiv.Comp,NumQuad) 
  fuDG  = zeros(FeDG.DoF,FeDG.DoF,NumQuad) 
  fuT  = zeros(FeT.DoF,FeT.Comp,NumQuad) 
  fuDivHDiv = zeros(FeHDiv.DoF,NumQuad) 
  fuGradDG = zeros(FeDG.DoF,FeDG.DoF,2,NumQuad) 

  #computation of the ansatz functions in the quadrature points
  for iQ = 1 : NumQuad
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fuT[iD,iComp,iQ] = FeT.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeHDiv.DoF
        fuHDiv[iD,iComp,iQ] = FeHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iD = 1 : FeDG.DoF
      fuGradDG[iD,1,1,iQ] = FeDG.Gradphi[iD,1,1](Points[iQ,1],Points[iQ,2])
      fuGradDG[iD,1,2,iQ] = FeDG.Gradphi[iD,1,2](Points[iQ,1],Points[iQ,2])
    end
    for iD = 1 : FeDG.DoF
      fuDG[iD,1,iQ] = FeDG.phi[iD,1](Points[iQ,1],Points[iQ,2])
    end
    for iD = 1 : FeHDiv.DoF
      fuDivHDiv[iD,iQ] = FeHDiv.Divphi[iD,1](Points[iQ,1],Points[iQ,2])
    end
  end

  #local variables
  cDGLoc = zeros(FeDG.DoF)
  uHDivLoc = zeros(FeHDiv.DoF)
  RhsLoc = zeros(FeT.DoF)

  uHDivLocLeft = zeros(FeHDiv.DoF)
  uHDivLocRight = zeros(FeHDiv.DoF)
  cDGLocLeft = zeros(FeDG.DoF)
  cDGLocRight = zeros(FeDG.DoF)

  uuGradDGLoc = zeros(2)
  uuHDivLoc = zeros(2)
 
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
 
  #copy from global variables to local variables  
  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iD = 1 : FeDG.DoF
      ind = FeDG.Glob[iD,iF]  
      cDGLoc[iD] = cDG[ind]
    end 
    @inbounds for iD = 1 : FeHDiv.DoF
      ind = FeHDiv.Glob[iD,iF]  
      uHDivLoc[iD] = uHDiv[ind]
    end  
    @. RhsLoc = 0.0

    for iQ = 1 : NumQuad
      #computation of local variables in a quadrature point
      uuDGLoc = 0.0
      @. uuGradDGLoc = 0.0
      @inbounds for iD = 1 : FeDG.DoF
        uuDGLoc += cDGLoc[iD] * fuDG[iD,1,iQ] 
        @views @. uuGradDGLoc += cDGLoc[iD] * fuGradDG[iD,1,:,iQ] 
      end
      uuDivHDivLoc = 0.0
      @. uuHDivLoc = 0.0
      @inbounds for iD = 1 : FeHDiv.DoF
        uuDivHDivLoc += uHDivLoc[iD] * fuDivHDiv[iD,iQ] 
        uuHDivLoc[1] += uHDivLoc[iD] * fuHDiv[iD,1,iQ]
        uuHDivLoc[2] += uHDivLoc[iD] * fuHDiv[iD,2,iQ] 
      end
      GradDGHDiv = uuGradDGLoc' * uuHDivLoc
      DGDivHDiv = uuDivHDivLoc * uuDGLoc
      innersum = Grid.Faces[iF].Orientation * (DGDivHDiv + GradDGHDiv)
      #product incoming functions and test function
      @inbounds for iD = 1 : FeT.DoF
        product = innersum * fuT[iD,1,iQ]
        RhsLoc[iD] += Weights[iQ] * product
      end 
    end

    #save the new without upwind 
    @inbounds for iD = 1 : FeT.DoF
      ind = FeT.Glob[iD,iF] 
      Rhs[ind] += RhsLoc[iD]
    end  
  end 

  #UPWIND on Edges
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd) 
  uFFRef  = zeros(FeHDiv.DoF, FeHDiv.Comp,NumQuadL,3)
  hFFRef  = zeros(FeDG.DoF,FeDG.DoF,NumQuadL,3)
  fTRef  = zeros(FeT.DoF,FeT.Comp,NumQuadL,3)
  uFFLocL = zeros(FeHDiv.DoF)
  uFFLocR = zeros(FeHDiv.DoF)

  #testfunctions DG scalar 
  for iQ = 1 : NumQuadL 
    for iComp = 1 : FeT.Comp
      for iD = 1 : FeT.DoF
        fTRef[iD,iComp,iQ,1] = FeT.phi[iD,iComp](PointsL[iQ],-1.0)
        fTRef[iD,iComp,iQ,2] = FeT.phi[iD,iComp](-PointsL[iQ],PointsL[iQ])
        fTRef[iD,iComp,iQ,3] = FeT.phi[iD,iComp](-1.0,-PointsL[iQ])
      end
    end
    #uf RT UHDiv
    for iComp = 1 : FeHDiv.Comp
      for iD = 1 : FeHDiv.DoF
        uFFRef[iD,iComp,iQ,1] = FeHDiv.phi[iD,iComp](PointsL[iQ],-1.0)
        uFFRef[iD,iComp,iQ,2] = FeHDiv.phi[iD,iComp](-PointsL[iQ],PointsL[iQ])
        uFFRef[iD,iComp,iQ,3] = FeHDiv.phi[iD,iComp](-1.0,-PointsL[iQ])
      end
    end
    #hf
    for iComp = 1 : FeDG.Comp
      for iD = 1 : FeDG.DoF
        hFFRef[iD,iComp,iQ,1] = FeDG.phi[iD,iComp](PointsL[iQ],-1.0)
        hFFRef[iD,iComp,iQ,2] = FeDG.phi[iD,iComp](-PointsL[iQ],PointsL[iQ])
        hFFRef[iD,iComp,iQ,3] = FeDG.phi[iD,iComp](-1.0,-PointsL[iQ])
      end
    end
  end

  #allocation Mloc
  MLoc11 = zeros(FeT.DoF, FeHDiv.DoF)
  MLoc12 = zeros(FeT.DoF, FeHDiv.DoF)
  MLoc21 = zeros(FeT.DoF, FeHDiv.DoF)
  MLoc22 = zeros(FeT.DoF, FeHDiv.DoF)
  
  Div = similar(Rhs)
  @. Div=0
  @. Rhs=0
  @inbounds for iE in 1:Grid.NumEdges
    Edge = Grid.Edges[iE]
    if length(Edge.F) > 1
      iFL = Edge.F[1]
      EdgeTypeL = Edge.FE[1]
      iFR = Edge.F[2]
      EdgeTypeR = Edge.FE[2]
      #computation normales of edges
      @views nBarLocL = Grid.nBar[:, EdgeTypeL]*Grid.Faces[iFL].Orientation
      @views nBarLocR = Grid.nBar[:, EdgeTypeR]*Grid.Faces[iFR].Orientation
      #gamma upwind value
      uE = 0.0 
      uE1 = 0.0
      @inbounds for iD = 1 : FeHDiv.DoF
        ind = FeHDiv.Glob[iD,iFL]  
        uHDivLocLeft[iD] = uHDiv[ind]
        ind = FeHDiv.Glob[iD,iFR]  
        uHDivLocRight[iD] = uHDiv[ind]
      end 
      @inbounds for iD = 1 : FeDG.DoF
        ind = FeDG.Glob[iD,iFL]  
        cDGLocLeft[iD] = cDG[ind]
        ind = FeDG.Glob[iD,iFR]  
        cDGLocRight[iD] = cDG[ind]
      end 
      @inbounds for iQ = 1:NumQuadL
        @inbounds for iD = 1 : FeHDiv.DoF
          @views uE += WeightsL[iQ]*(uHDivLocLeft[iD]*(uFFRef[iD,:,iQ,EdgeTypeL]'*nBarLocL))
          @views uE1 += WeightsL[iQ]*(uHDivLocRight[iD]*(uFFRef[iD,:,iQ,EdgeTypeR]'*nBarLocR))
        end
      end
      #upwind value
      gammaU = 0.5
      gammaLoc = uE > 0 ? gammaU : -gammaU

      @. MLoc11 = 0
      @. MLoc12 = 0
      @. MLoc21 = 0
      @. MLoc22 = 0
      
      @inbounds for iQ in 1:NumQuadL
        hFFL = 0.0
        hFFR = 0.0
        @inbounds for iD in 1:FeHDiv.DoF
          @views uFFLocL[iD] = uFFRef[iD,:,iQ,EdgeTypeL]' * nBarLocL
          @views uFFLocR[iD] = uFFRef[iD,:,iQ,EdgeTypeR]' * nBarLocR
        end
        @inbounds for iD in 1:FeDG.DoF
          hFFL += hFFRef[iD,1,iQ,EdgeTypeL] * cDGLocLeft[iD]
          hFFR += hFFRef[iD,1,iQ,EdgeTypeR] * cDGLocRight[iD]
        end

        @inbounds for iDoF  = 1 : FeT.DoF 
          @inbounds for jDoF  = 1 : FeHDiv.DoF
            MLoc11[iDoF,jDoF] += WeightsL[iQ] * (fTRef[iDoF,1,iQ,EdgeTypeL] * hFFL) * uFFLocL[jDoF]'
            MLoc12[iDoF,jDoF] += WeightsL[iQ] * (fTRef[iDoF,1,iQ,EdgeTypeL] * hFFR) * uFFLocR[jDoF]'
            MLoc21[iDoF,jDoF] += WeightsL[iQ] * (fTRef[iDoF,1,iQ,EdgeTypeR] * hFFL) * uFFLocL[jDoF]'
            MLoc22[iDoF,jDoF] += WeightsL[iQ] * (fTRef[iDoF,1,iQ,EdgeTypeR] * hFFR) * uFFLocR[jDoF]'
          end
        end
      end
      #all together over edges
      @inbounds for iD = 1 : FeT.DoF
        indL = FeT.Glob[iD,iFL] 
        @views Rhs[indL] += (+0.5 - gammaLoc) * MLoc11[iD,:]' * uHDivLocLeft
        @views Rhs[indL] += (-0.5 + gammaLoc) * MLoc12[iD,:]' * uHDivLocRight
        indR = FeT.Glob[iD,iFR] 
        @views Rhs[indR] += (+0.5 + gammaLoc) * MLoc21[iD,:]' * uHDivLocLeft
        @views Rhs[indR] += (-0.5 - gammaLoc) * MLoc22[iD,:]' * uHDivLocRight

        Div[indL] += uHDiv[iE]  * Grid.Faces[iFL].OrientE[EdgeTypeL] 
        Div[indR] += uHDiv[iE]  * Grid.Faces[iFR].OrientE[EdgeTypeR]
        if indL == 2 || indR == 2
        @show Rhs[2], Div[2]
        @show uE, uE1
        @show MLoc11
      
       
        end
      end
    end
  end
  #stop
  #@. Rhs=Div
end
