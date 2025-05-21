function DivMomentumVector!(backend,FTB,Rhs,FeTHDiv::HDivElement,uHDiv,FeHDiv::HDivElement,
  uVecDG,FeVecDG::VectorElement,Grid,ElemType::Grids.ElementType,QuadOrd,Jacobi)

  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  fuHDiv  = zeros(FeHDiv.DoF,FeHDiv.Comp,NumQuad)
  fuVecDG  = zeros(FeVecDG.DoF,FeVecDG.Comp,NumQuad)
  fuTHDiv  = zeros(FeTHDiv.DoF,FeTHDiv.Comp,NumQuad)
  fuDivHDiv = zeros(FeHDiv.DoF,NumQuad)
  fuGradVecDG = zeros(FeVecDG.DoF,FeVecDG.Comp,2,NumQuad)

  #computation of the ansatz functions in the quadrature points
  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : FeTHDiv.Comp
      @inbounds for iD = 1 : FeTHDiv.DoF
        fuTHDiv[iD,iComp,iQ] = FeTHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : FeTHDiv.Comp
      @inbounds for iD = 1 : FeHDiv.DoF
        fuHDiv[iD,iComp,iQ] = FeHDiv.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iD = 1 : FeVecDG.DoF
      fuVecDG[iD,1,iQ] = FeVecDG.phi[iD,1](Points[iQ,1],Points[iQ,2])
      fuVecDG[iD,2,iQ] = FeVecDG.phi[iD,2](Points[iQ,1],Points[iQ,2])
      fuVecDG[iD,3,iQ] = FeVecDG.phi[iD,3](Points[iQ,1],Points[iQ,2])
    end
    @inbounds for iD = 1 : FeHDiv.DoF
      fuDivHDiv[iD,iQ] = FeHDiv.Divphi[iD,1](Points[iQ,1],Points[iQ,2])
    end
    
    @inbounds for iD = 1 : FeVecDG.DoF
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
      @. innersum = (VecDGDivHDiv + GradVecDGHDiv)
      #computation of Jacobi
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      @inbounds for iD = 1 : FeTHDiv.DoF
        product = 0.0
        @inbounds for i = 1 : 3
          temp = 0.0  
          @inbounds for j = 1 : 2  
            temp +=  DF[i,j] * fuTHDiv[iD,j,iQ] / detDFLoc 
          end  
          product += innersum[i] * temp
        end  
        RhsLoc[iD] += product * Weights[iQ]
      end 
    end
    @inbounds for iD = 1 : FeTHDiv.DoF
      ind = FeTHDiv.Glob[iD,iF]
      Rhs[ind] -= RhsLoc[iD]
    end  
  end 

  #UPWIND on Edges
  NumQuadL, WeightsL, PointsL = QuadRule(Grids.Line(),QuadOrd)
  #uFFLocL = zeros(FeHDiv.DoF)
  #uFFLocR = zeros(FeHDiv.DoF)
  uHDivLocLeft = zeros(FeHDiv.DoF)
  uHDivLocRight = zeros(FeHDiv.DoF)
  uVecDGLocLeft = zeros(FeVecDG.DoF)
  uVecDGLocRight = zeros(FeVecDG.DoF)

  #distinction between Tri or Quad 
  if ElemType == Grids.Tri()
    PointsE = zeros(2,NumQuadL,3)
    @inbounds for iQ = 1 : NumQuadL
      PointsE[1,iQ,1] = PointsL[iQ]
      PointsE[2,iQ,1] = -1.0
      PointsE[1,iQ,2] = -PointsL[iQ]
      PointsE[2,iQ,2] = PointsL[iQ]
      PointsE[1,iQ,3] = -1.0
      PointsE[2,iQ,3] = PointsL[iQ]
    end
  elseif ElemType == Grids.Quad()
    PointsE = zeros(2,NumQuadL,4)
    @inbounds for iQ = 1 : NumQuadL
      PointsE[1,iQ,1] = PointsL[iQ]
      PointsE[2,iQ,1] = -1.0
      PointsE[1,iQ,2] = 1.0
      PointsE[2,iQ,2] = PointsL[iQ]
      PointsE[1,iQ,3] = PointsL[iQ]
      PointsE[2,iQ,3] = 1.0
      PointsE[1,iQ,4] = -1.0
      PointsE[2,iQ,4] = PointsL[iQ]
    end
  else 
    @show "wrong ElemType"
  end
  uFFRef  = zeros(FeHDiv.DoF, FeHDiv.Comp,NumQuadL,size(PointsE,3))
  uVecDGRef  = zeros(FeVecDG.DoF,FeVecDG.Comp,NumQuadL,size(PointsE,3))
  fTRef  = zeros(FeTHDiv.DoF,FeTHDiv.Comp,NumQuadL,size(PointsE,3))
  @inbounds for i = 1 : size(PointsE,3)
    @inbounds for iQ = 1 : NumQuadL 
      @inbounds for iComp = 1 : FeTHDiv.Comp
        @inbounds for iD = 1 : FeTHDiv.DoF
          fTRef[iD,iComp,iQ,i] = FeTHDiv.phi[iD,iComp](PointsE[1,iQ,i],PointsE[2,iQ,i])
        end
      end
      @inbounds for iComp = 1 : FeHDiv.Comp
        @inbounds for iD = 1 : FeHDiv.DoF
          uFFRef[iD,iComp,iQ,i] = FeHDiv.phi[iD,iComp](PointsE[1,iQ,i],PointsE[2,iQ,i])
        end
      end
      @inbounds for iComp = 1 : FeVecDG.Comp
        @inbounds for iD = 1 : FeVecDG.DoF
          uVecDGRef[iD,iComp,iQ,i] = FeVecDG.phi[iD,iComp](PointsE[1,iQ,i],PointsE[2,iQ,i])
        end
      end
    end
  end

  #allocation Mloc
  ILL = zeros(FeTHDiv.DoF)
  ILR = zeros(FeTHDiv.DoF)
  IRL = zeros(FeTHDiv.DoF)
  IRR = zeros(FeTHDiv.DoF)
  fTLocL = zeros(3,FeTHDiv.DoF)
  fTLocR = zeros(3,FeTHDiv.DoF)

  uL = zeros(3)
  uR = zeros(3)
  gammaU = 0.5
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
      uER = 0.0 
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
      @inbounds for iQ = 1:NumQuadL
        @inbounds for iD = 1 : FeHDiv.DoF
          @views uE += WeightsL[iQ]*(uHDivLocLeft[iD]*uFFRef[iD,:,iQ,EdgeTypeL]'*nBarLocL) 
        end
      end
      gammaLoc = uE > 0 ? gammaU : -gammaU

      @. ILL = 0
      @. ILR = 0
      @. IRL = 0
      @. IRR = 0

      @inbounds for iQ in 1: NumQuadL
        #computation of Jacobi EdgetypeL
        Jacobi(DFL,detDFL,pinvDFL,XL,Grid.Type,PointsE[1,iQ,EdgeTypeL],PointsE[2,iQ,EdgeTypeL],Grid.Faces[iFL], Grid)
        detDFLLoc = detDFL[1] * Grid.Faces[iFL].Orientation
        #EdgeTypeR
        Jacobi(DFR,detDFR,pinvDFR,XR,Grid.Type,PointsE[1,iQ,EdgeTypeR],PointsE[2,iQ,EdgeTypeR],Grid.Faces[iFR], Grid)
        detDFRLoc = detDFR[1] * Grid.Faces[iFR].Orientation

        @inbounds for iDoF  = 1 : FeTHDiv.DoF 
          @inbounds for j = 1 : 3
            fTLocL[j,iDoF] = (1/detDFLLoc) * (DFL[j,1] * fTRef[iDoF,1,iQ,EdgeTypeL] +
              DFL[j,2] * fTRef[iDoF,2,iQ,EdgeTypeL])
            fTLocR[j,iDoF] = (1/detDFRLoc) * (DFR[j,1] * fTRef[iDoF,1,iQ,EdgeTypeR] +
              DFR[j,2] * fTRef[iDoF,2,iQ,EdgeTypeR])
          end
        end
        uFFLocL = 0.0
        uFFLocR = 0.0
        @inbounds for iD in 1:FeHDiv.DoF
          @views uFFLocL = uHDivLocLeft[iD] * uFFRef[iD,:,iQ,EdgeTypeL]' * nBarLocL
          @views uFFLocR = uHDivLocRight[iD] * uFFRef[iD,:,iQ,EdgeTypeR]' * nBarLocR
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
              ILL[iDoF] += WeightsL[iQ] * (fTLocL[k,iDoF]' * uL[k])
              ILR[iDoF] += WeightsL[iQ] * (fTLocL[k,iDoF]' * uR[k])
              IRL[iDoF] += WeightsL[iQ] * (fTLocR[k,iDoF]' * uL[k])
              IRR[iDoF] += WeightsL[iQ] * (fTLocR[k,iDoF]' * uR[k])
            end
          end
          ILL[iDoF] *= uFFLocL
          ILR[iDoF] *= uFFLocL
          IRL[iDoF] *= uFFLocR
          IRR[iDoF] *= uFFLocR
        end
      end
      #all together over edges
      @inbounds for iD = 1 : FeTHDiv.DoF
#       indR = FeTHDiv.Glob[iD,iFR]
#       Rhs[indR] = 0.5 * (IRR[iD] - IRL[iD])
#       indL = FeTHDiv.Glob[iD,iFL]
#       Rhs[indL] = 0.5 * (ILL[iD] - ILR[iD])
        
        if gammaLoc > 0.0
          indR = FeTHDiv.Glob[iD,iFR]
          Rhs[indR] = -1.0 * IRL[iD]
          Rhs[indR] = +1.0 * IRR[iD]
        else
          indL = FeTHDiv.Glob[iD,iFL]
          Rhs[indL] = +1.0  * ILL[iD]
          Rhs[indL] = -1.0  * ILR[iD]
        end  
      end
    end
  end
end
