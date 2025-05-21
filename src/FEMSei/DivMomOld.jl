function DivMomentumVectorOld!(backend,FTB,Rhs,FeTHDiv::HDivElement,uHDiv,FeHDiv::HDivElement,
  uVecDG,FeVecDG::VectorElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  fuHDiv  = zeros(FeHDiv.DoF,FeHDiv.Comp,NumQuad)
  fuVecDG  = zeros(FeVecDG.DoF,FeVecDG.DoF,NumQuad)
  fuTHDiv  = zeros(FeTHDiv.DoF,FeTHDiv.Comp,NumQuad)
  fuDivHDiv = zeros(FeHDiv.DoF,NumQuad)
  fuGradVecDG = zeros(FeVecDG.DoF,FeVecDG.DoF,2,NumQuad)

  #computation of the ansatz functions in the quadrature points
  for iQ = 1 : NumQuad
    for iComp = 1 : FeTHDiv.Comp
      for iD = 1 : FeTHDiv.DoF
        fuTHDiv[iD,iComp,iQ] = FeTHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iComp = 1 : FeTHDiv.Comp
      for iD = 1 : FeHDiv.DoF
        fuHDiv[iD,iComp,iQ] = FeHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    for iD = 1 : FeVecDG.DoF
      fuVecDG[iD,1,iQ] = FeVecDG.phi[iD,1](Points[iQ,1],Points[iQ,2])
      fuVecDG[iD,2,iQ] = FeVecDG.phi[iD,2](Points[iQ,1],Points[iQ,2])
      fuVecDG[iD,3,iQ] = FeVecDG.phi[iD,3](Points[iQ,1],Points[iQ,2])
    end
    for iD = 1 : FeHDiv.DoF
      fuDivHDiv[iD,iQ] = FeHDiv.Divphi[iD,1](Points[iQ,1],Points[iQ,2])
    end
    
    for iD = 1 : FeVecDG.DoF
      fuGradVecDG[iD,1,1,iQ] = FeVecDG.Gradphi[iD,1,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,1,2,iQ] = FeVecDG.Gradphi[iD,1,2](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,2,1,iQ] = FeVecDG.Gradphi[iD,2,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,2,2,iQ] = FeVecDG.Gradphi[iD,2,2](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,3,1,iQ] = FeVecDG.Gradphi[iD,3,1](Points[iQ,1],Points[iQ,2])
      fuGradVecDG[iD,3,2,iQ] = FeVecDG.Gradphi[iD,3,2](Points[iQ,1],Points[iQ,2])
    end
  end


  #local variables
  uVecDGLoc = zeros(FeVecDG.DoF) 
  uuVecDGLoc = zeros(3) 
  uuGradVecDGLoc = zeros(3,2) 
  uHDivLoc = zeros(FeHDiv.DoF)
  uuHDivLoc = zeros(2)
  RhsLoc = zeros(FeTHDiv.DoF)

  GradVecDGHDiv = zeros(3)
  VecDGDivHDiv = zeros(3)
  
  #Jacobi-Allocation
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  DFL = zeros(3,2)
  detDFL = zeros(1)
  pinvDFL = zeros(3,2)
  XL = zeros(3)

  DFR = zeros(3,2)
  detDFR = zeros(1)
  pinvDFR = zeros(3,2)
  XR = zeros(3)
  innersum = zeros(3)

  #copy from global variables to local variables  
  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iD = 1 : FeVecDG.DoF
      ind = FeVecDG.Glob[iD,iF]  
      uVecDGLoc[iD] = uVecDG[ind]
    end  
    @inbounds for iD = 1 : FeHDiv.DoF
      ind = FeHDiv.Glob[iD,iF]  
      uHDivLoc[iD] = uHDiv[ind]
    end  
    @. RhsLoc = 0.0

    @inbounds for iQ = 1 : NumQuad
      @. uuVecDGLoc = 0.0
      @. uuGradVecDGLoc = 0.0
      #computation of local variables in a quadrature point
      @inbounds for iD = 1 : FeVecDG.DoF
        @inbounds for i = 1 : 3
          uuVecDGLoc[i] += uVecDGLoc[iD] * fuVecDG[iD,i,iQ] 
          @inbounds for j = 1 : 2  
            uuGradVecDGLoc[i,j] += uVecDGLoc[iD] * fuGradVecDG[iD,i,j,iQ]  
          end
        end
      end
      uuDivHDivLoc = 0.0
      @. uuHDivLoc = 0.0
      @inbounds for iD = 1 : FeHDiv.DoF
        uuDivHDivLoc += uHDivLoc[iD] * fuDivHDiv[iD,iQ] 
        @inbounds for i = 1 : 2
          uuHDivLoc[i] += uHDivLoc[iD] * fuHDiv[iD,i,iQ]
        end 
      end

      @inbounds for i = 1 : 3
        GradVecDGHDiv[i] = 0.0  
        VecDGDivHDiv[i] = uuDivHDivLoc * uuVecDGLoc[i]
        @inbounds for j = 1 : 2  
          GradVecDGHDiv[i] += uuGradVecDGLoc[i,j] * uuHDivLoc[j]
        end
      end  
      @. innersum = Grid.Faces[iF].Orientation * (VecDGDivHDiv + GradVecDGHDiv)
      #computation of Jacobi
      Jacobi!(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      product = 0.0
      @inbounds for iD = 1 : FeTHDiv.DoF
        @inbounds for i = 1 : 3
          temp = 0.0  
          @inbounds for j = 1 : 2  
            temp +=  DF[i,j] * fuTHDiv[iD,j,iQ]  
          end  
          product += innersum[i] * temp
        end  
      end 
    end
    @inbounds for iD = 1 : FeTHDiv.DoF
      ind = FeTHDiv.Glob[iD,iF]
      Rhs[ind] += RhsLoc[iD]
    end  
  end 

  #UPWIND on Edges
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  uFFRef  = zeros(FeHDiv.DoF, FeHDiv.Comp,NumQuadL,4)
  uVecDGRef  = zeros(FeVecDG.DoF,FeVecDG.Comp,NumQuadL,4)
  fTRef  = zeros(FeTHDiv.DoF,FeTHDiv.Comp,NumQuadL,4)
  uFFLocL = zeros(FeHDiv.DoF)
  uFFLocR = zeros(FeHDiv.DoF)
  uHDivLocLeft = zeros(FeHDiv.DoF)
  uHDivLocRight = zeros(FeHDiv.DoF)
  uVecDGLocLeft = zeros(FeVecDG.DoF)
  uVecDGLocRight = zeros(FeVecDG.DoF)
  
  for iQ = 1 : NumQuadL 
    for iComp = 1 : FeTHDiv.Comp
      for iD = 1 : FeTHDiv.DoF
        fTRef[iD,iComp,iQ,1] = FeTHDiv.phi[iD,iComp](PointsL[iQ],-1.0)
        fTRef[iD,iComp,iQ,2] = FeTHDiv.phi[iD,iComp](1.0,PointsL[iQ])
        fTRef[iD,iComp,iQ,3] = FeTHDiv.phi[iD,iComp](PointsL[iQ],1.0)
        fTRef[iD,iComp,iQ,4] = FeTHDiv.phi[iD,iComp](-1.0,PointsL[iQ])
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
    for iComp = 1 : FeVecDG.Comp
      for iD = 1 : FeVecDG.DoF
        uVecDGRef[iD,iComp,iQ,1] = FeVecDG.phi[iD,iComp](PointsL[iQ],-1.0)
        uVecDGRef[iD,iComp,iQ,2] = FeVecDG.phi[iD,iComp](1.0,PointsL[iQ])
        uVecDGRef[iD,iComp,iQ,3] = FeVecDG.phi[iD,iComp](PointsL[iQ],1.0)
        uVecDGRef[iD,iComp,iQ,4] = FeVecDG.phi[iD,iComp](-1.0,PointsL[iQ])
      end
    end
  end

  #allocation Mloc
  MLoc11 = zeros(FeTHDiv.DoF, FeHDiv.DoF)
  MLoc12 = zeros(FeTHDiv.DoF, FeHDiv.DoF)
  MLoc21 = zeros(FeTHDiv.DoF, FeHDiv.DoF)
  MLoc22 = zeros(FeTHDiv.DoF, FeHDiv.DoF)
  fTLocL = zeros(3,FeTHDiv.DoF)
  fTLocR = zeros(3,FeTHDiv.DoF)

  uL = zeros(3)
  uR = zeros(3)

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
      @inbounds for iD = 1 : FeVecDG.DoF
        ind = FeVecDG.Glob[iD,iFL]  
        uVecDGLocLeft[iD] = uVecDG[ind]
        ind = FeVecDG.Glob[iD,iFR]  
        uVecDGLocRight[iD] = uVecDG[ind]
      end 
      @inbounds for iQ = 1:length(WeightsL)
        @inbounds for iD = 1 : FeHDiv.DoF
          @views uE += WeightsL[iQ]*(uHDivLocLeft[iD]*(uFFRef[iD,:,iQ,EdgeTypeL]'*nBarLocL) +
            uHDivLocRight[iD]*(uFFRef[iD,:,iQ,EdgeTypeR]'*nBarLocR))
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
        #computation of Jacobi EdgetypeL
        Jacobi!(DFL,detDFL,pinvDFL,XL,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iFL], Grid)
        detDFLLoc = detDFL[1] 
        #EdgeTypeR
        Jacobi!(DFR,detDFR,pinvDFR,XR,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iFR], Grid)
        detDFRLoc = detDFR[1]

        @inbounds for iDoF  = 1 : FeTHDiv.DoF 
          @inbounds for j = 1 : 3
            fTLocL[j,iDoF] = (1/detDFLLoc) * (DFL[j,1] * fTRef[iDoF,1,iQ,EdgeTypeL] +
              DFL[j,2] * fTRef[iDoF,2,iQ,EdgeTypeL])
            fTLocR[j,iDoF] = (1/detDFRLoc) * (DFR[j,1] * fTRef[iDoF,1,iQ,EdgeTypeR] +
              DFR[j,2] * fTRef[iDoF,2,iQ,EdgeTypeR])
          end
        end
        @inbounds for iD in 1:FeHDiv.DoF
          @views uFFLocL[iD] = uFFRef[iD,:,iQ,EdgeTypeL]' * nBarLocL
          @views uFFLocR[iD] = uFFRef[iD,:,iQ,EdgeTypeR]' * nBarLocR
        end
        @. uL = 0.0
        @. uR = 0.0
      
        @inbounds for iD in 1:FeVecDG.DoF 
          uL[1] +=  uVecDGRef[iD,1,iQ,EdgeTypeL] * uVecDGLocLeft[iD]          
          uL[2] +=  uVecDGRef[iD,2,iQ,EdgeTypeL] * uVecDGLocLeft[iD]          
          uL[3] +=  uVecDGRef[iD,3,iQ,EdgeTypeL] * uVecDGLocLeft[iD]          
          uR[1] +=  uVecDGRef[iD,1,iQ,EdgeTypeR] * uVecDGLocRight[iD]
          uR[2] +=  uVecDGRef[iD,2,iQ,EdgeTypeR] * uVecDGLocRight[iD]
          uR[3] +=  uVecDGRef[iD,3,iQ,EdgeTypeR] * uVecDGLocRight[iD]
        end
        @inbounds for iDoF  = 1 : FeTHDiv.DoF 
          @inbounds for jDoF  = 1 : FeHDiv.DoF 
            @inbounds for k = 1 : 3  
              MLoc11[iDoF,jDoF] += WeightsL[iQ] * (fTLocL[k,iDoF]' * uL[k]) * uFFLocL[jDoF]
              MLoc12[iDoF,jDoF] += WeightsL[iQ] * (fTLocL[k,iDoF]' * uR[k]) * uFFLocR[jDoF]
              MLoc21[iDoF,jDoF] += WeightsL[iQ] * (fTLocR[k,iDoF]' * uL[k]) * uFFLocL[jDoF]
              MLoc22[iDoF,jDoF] += WeightsL[iQ] * (fTLocR[k,iDoF]' * uR[k]) * uFFLocR[jDoF]
            end
          end
        end
      end
      #all together over edges
      @inbounds for iD = 1 : FeTHDiv.DoF
        indL = FeTHDiv.Glob[iD,iFL] 
        @views Rhs[indL] += (-0.5 - gammaLoc) * MLoc11[iD,:]' * uHDivLocLeft
        @views Rhs[indL] += (-0.5 + gammaLoc) * MLoc12[iD,:]' * uHDivLocRight
        indR = FeTHDiv.Glob[iD,iFR] 
        @views Rhs[indR] += (+0.5 + gammaLoc) * MLoc21[iD,:]' * uHDivLocLeft
        @views Rhs[indR] += (+0.5 - gammaLoc) * MLoc22[iD,:]' * uHDivLocRight
      end
    end
  end
end
