function TopographySmoothing2!(hFCG,hCG,CG,Global,HyperDDiv)

  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  Div = zeros(CG.NumG)
  hF = zeros(CG.NumG)
  DivCG= zeros(OP,OP)
  @. Div = 0.0
  @. hF = 0.0
  D1cCG = zeros(OP,OP)
  D2cCG = zeros(OP,OP) 
  grad1CG = zeros(OP,OP) 
  grad2CG = zeros(OP,OP) 

  D1gradCG = D1cCG
  D2gradCG = D2cCG

  vC1 = grad1CG
  vC2 = grad2CG

  # Hyperdiffusion 
  @inbounds for iF = 1:NF
    @views J = Global.Metric.J[:,:,:,1,iF];
    @views dXdxI = Global.Metric.dXdxI[:,:,:,1,:,:,iF]

    @views mul!(D1cCG[:,:],CG.DS,hCG[:,:,iF])
    @views mul!(D2cCG[:,:],hCG[:,:,iF],CG.DST)

    @views @. grad1CG[:,:] = (dXdxI[:,:,1,1,1] + dXdxI[:,:,2,1,1]) * D1cCG[:,:] + 
      (dXdxI[:,:,1,2,1] + dXdxI[:,:,2,2,1]) * D2cCG[:,:]
    @views @. grad2CG[:,:] = (dXdxI[:,:,1,1,2] + dXdxI[:,:,2,1,2]) * D1cCG[:,:] + 
      (dXdxI[:,:,1,2,2] + dXdxI[:,:,2,2,2]) * D2cCG[:,:]

    @views @. D1gradCG[:,:] = (dXdxI[:,:,1,1,1] + dXdxI[:,:,2,1,1]) * grad1CG[:,:] + 
      (dXdxI[:,:,1,1,2] + dXdxI[:,:,2,1,2]) * grad2CG[:,:]
    @views @. D2gradCG[:,:] = (dXdxI[:,:,1,1,2,1] + dXdxI[:,:,2,1,2,1]) * grad1CG[:,:] + 
      (dXdxI[:,:,1,2,2] + dXdxI[:,:,2,2,2]) * grad2CG[:,:]

    @views mul!(vC1[:,:],CG.DW,D1gradCG[:,:])
    @views mul!(vC2[:,:],D2gradCG[:,:],CG.DWT)
    @views @. DivCG[:,:] = (vC1[:,:] + vC2[:,:]) / (J[:,:,1] + J[:,:,2])
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        Div[ind] += DivCG[iP,jP] / CG.M[1,ind]
      end
    end
  end

  # Hyperdiffusion 
  @inbounds for iF = 1:NF
    @views J = Global.Metric.J[:,:,:,1,iF];
    @views dXdxI = Global.Metric.dXdxI[:,:,:,:,:,:,iF]
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        DivCG[iP,jP] = Div[ind]
      end
    end
    @views mul!(D1cCG[:,:],CG.DS,DivCG[:,:])
    @views mul!(D2cCG[:,:],DivCG[:,:],CG.DST)
  
    @views @. grad1CG[:,:] = (dXdxI[:,:,1,1,1] + dXdxI[:,:,2,1,1]) * D1cCG[:,:] + 
      (dXdxI[:,:,1,2,1] + dXdxI[:,:,2,2,1]) * D2cCG[:,:]
    @views @. grad2CG[:,:] = (dXdxI[:,:,1,1,2] + dXdxI[:,:,2,1,2]) * D1cCG[:,:] + 
      (dXdxI[:,:,1,2,2] + dXdxI[:,:,2,2,2]) * D2cCG[:,:]

    @views @. D1gradCG[:,:] = (dXdxI[:,:,1,1,1] + dXdxI[:,:,2,1,1]) * grad1CG[:,:] + 
      (dXdxI[:,:,1,1,2] + dXdxI[:,:,2,1,2]) * grad2CG[:,:]
    @views @. D2gradCG[:,:] = (dXdxI[:,:,1,1,2,1] + dXdxI[:,:,2,1,2,1]) * grad1CG[:,:] + 
      (dXdxI[:,:,1,2,2] + dXdxI[:,:,2,2,2]) * grad2CG[:,:]

    @views mul!(vC1[:,:],CG.DW,D1gradCG[:,:])
    @views mul!(vC2[:,:],D2gradCG[:,:],CG.DWT)
    @views @. DivCG[:,:] -= HyperDDiv*(vC1[:,:] + vC2[:,:]) / (J[:,:,1] + J[:,:,2])

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        hF[ind] += DivCG[iP,jP] / CG.M[1,ind]
      end
    end
  end  
  ExchangeData!(hF,Global.Exchange)

  @inbounds for iF = 1:NF
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        hFCG[iP,jP,iF] = hF[ind]
      end
    end
  end
end

function TopographySmoothing1!(hFCG,hCG,CG,Global,HyperDDiv)

  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  J = Global.Metric.J
  Div = zeros(CG.NumG)
  hF = zeros(CG.NumG)
  DivCG= zeros(OP,OP)
  @. Div = 0.0
  @. hF = 0.0
  D1cCG = zeros(OP,OP)
  D2cCG = zeros(OP,OP) 
  grad1CG = zeros(OP,OP) 
  grad2CG = zeros(OP,OP) 

  D1gradCG = D1cCG
  D2gradCG = D2cCG

  vC1 = grad1CG
  vC2 = grad2CG

  # Diffusion 
  @inbounds for iF = 1:NF
    @views J = Global.Metric.J[:,:,:,1,iF];
    @views dXdxI = Global.Metric.dXdxI[:,:,:,1,:,:,iF]

    @views mul!(D1cCG[:,:],CG.DS,hCG[:,:,iF])
    @views mul!(D2cCG[:,:],hCG[:,:,iF],CG.DST)

  
    @views @. grad1CG[:,:] = (dXdxI[:,:,1,1,1] + dXdxI[:,:,2,1,1]) * D1cCG[:,:] + 
      (dXdxI[:,:,1,2,1] + dXdxI[:,:,2,2,1]) * D2cCG[:,:]
    @views @. grad2CG[:,:] = (dXdxI[:,:,1,1,2] + dXdxI[:,:,2,1,2]) * D1cCG[:,:] + 
      (dXdxI[:,:,1,2,2] + dXdxI[:,:,2,2,2]) * D2cCG[:,:]

    @views @. D1gradCG[:,:] = (dXdxI[:,:,1,1,1] + dXdxI[:,:,2,1,1]) * grad1CG[:,:] + 
      (dXdxI[:,:,1,1,2] + dXdxI[:,:,2,1,2]) * grad2CG[:,:]
    @views @. D2gradCG[:,:] = (dXdxI[:,:,1,2,1] + dXdxI[:,:,2,2,1]) * grad1CG[:,:] + 
      (dXdxI[:,:,1,2,2] + dXdxI[:,:,2,2,2]) * grad2CG[:,:]

    @views mul!(vC1[:,:],CG.DW,D1gradCG[:,:])
    @views mul!(vC2[:,:],D2gradCG[:,:],CG.DWT)
    @views @. DivCG[:,:] = HyperDDiv*(vC1[:,:] + vC2[:,:]) / (J[:,:,1] + J[:,:,2])
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        hF[ind] += DivCG[iP,jP] / CG.M[1,ind]
      end
    end
  end
  ExchangeData!(hF,Global.Exchange)

  @inbounds for iF = 1:NF
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        hFCG[iP,jP,iF] = hF[ind]
      end
    end
  end
end
