function CrossRhs!(backend,FTB,Cross,FeT::HDivElement,q,qFeF::ScalarElement,u,uFeF::HDivElement,Grid,
  ElemType::Grids.ElementType,QuadOrd,Jacobi)
  NumQuad, Weights, Points = QuadRule(ElemType,QuadOrd)
  qFFRef  = zeros(qFeF.Comp,qFeF.DoF,NumQuad)
  fTRef  = zeros(FeT.Comp,FeT.DoF,NumQuad)

  @inbounds for iQ = 1 : NumQuad
    @inbounds for iComp = 1 : FeT.Comp
      @inbounds for iD = 1 : FeT.DoF
        fTRef[iComp,iD,iQ] = FeT.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
    @inbounds for iComp = 1 : qFeF.Comp
      @inbounds for iD = 1 : qFeF.DoF
        qFFRef[iComp,iD,iQ] = qFeF.phi[iD,iComp](Points[iQ,1],Points[iQ,2])
      end
    end
  end
  uFFRef = fTRef

  CrossLoc = zeros(FeT.DoF)
  uLoc = zeros(uFeF.DoF)
  qLoc = zeros(qFeF.DoF)
  cu1 = zeros(NumQuad)
  cu2 = zeros(NumQuad)
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)

  @inbounds for iF = 1 : Grid.NumFaces
    @. CrossLoc = 0
    @inbounds for iDoF = 1 : uFeF.DoF
      ind = uFeF.Glob[iDoF,iF]  
      uLoc[iDoF] = u[ind]
    end  
    @inbounds for iDoF = 1 : qFeF.DoF
      ind = qFeF.Glob[iDoF,iF]  
      qLoc[iDoF] = q[ind]
    end  
    @inbounds for iQ = 1 : NumQuad
      Jacobi(DF,detDF,pinvDF,X,Grid.Type,Points[iQ,1],Points[iQ,2],Grid.Faces[iF], Grid)
      detDFLoc = detDF[1]
      uFLoc1 = 0.0
      uFLoc2 = 0.0
      @inbounds for iDoF = 1 : uFeF.DoF
        uFLoc1 += uFFRef[1,iDoF,iQ] * uLoc[iDoF]
        uFLoc2 += uFFRef[2,iDoF,iQ] * uLoc[iDoF]
      end
      u1 = (DF[1,1] * uFLoc1 + DF[1,2] * uFLoc2) 
      u2 = (DF[2,1] * uFLoc1 + DF[2,2] * uFLoc2)
      u3 = (DF[3,1] * uFLoc1 + DF[3,2] * uFLoc2)
      
      qFLoc,k= Coriolis(X,Grid.Form)

      kcru1 = k[2] * u3 - k[3] * u2
      kcru2 = -(k[1] * u3 - k[3] * u1)
      kcru3 = k[1] * u2 - k[2] * u1
      cu1 = (DF[1,1] * kcru1 + DF[2,1] * kcru2 + DF[3,1] * kcru3)
      cu2 = (DF[1,2] * kcru1 + DF[2,2] * kcru2 + DF[3,2] * kcru3)    
      @inbounds for iDoF = 1 : qFeF.DoF
        qFLoc += qFFRef[1,iDoF,iQ] * qLoc[iDoF]
      end
      @inbounds for iDoF = 1 : FeT.DoF
        CrossLoc[iDoF] -= Weights[iQ] * qFLoc * (fTRef[1,iDoF,iQ] * cu1 +
          fTRef[2,iDoF,iQ] * cu2) / detDFLoc
      end    
    end  
    @inbounds for iDoF = 1 : FeT.DoF
      ind = FeT.Glob[iDoF,iF]
      Cross[ind] += CrossLoc[iDoF]
    end
  end
end
