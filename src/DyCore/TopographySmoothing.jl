function TopographySmoothing2!(hFCG,hCG,CG,Global,HyperDDiv)

  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  JF = Global.Metric.JF
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
    @views JC = Global.Metric.JC[:,:,1,iF];
    @views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]

    @views mul!(D1cCG[:,:],CG.DS,hCG[:,:,iF])
    @views mul!(D2cCG[:,:],hCG[:,:,iF],CG.DST)

    @views @. grad1CG[:,:] = dXdxIC[:,:,1,1,1] * D1cCG[:,:] + dXdxIC[:,:,1,2,1] * D2cCG[:,:]
    @views @. grad2CG[:,:] = dXdxIC[:,:,1,1,2] * D1cCG[:,:] + dXdxIC[:,:,1,2,2] * D2cCG[:,:]

    @views @. D1gradCG[:,:] = dXdxIC[:,:,1,1,1] * grad1CG[:,:] + dXdxIC[:,:,1,1,2] * grad2CG[:,:]
    @views @. D2gradCG[:,:] = dXdxIC[:,:,1,2,1] * grad1CG[:,:] + dXdxIC[:,:,1,2,2] * grad2CG[:,:]

    @views mul!(vC1[:,:],CG.DW,D1gradCG[:,:])
    @views mul!(vC2[:,:],D2gradCG[:,:],CG.DWT)
    @views @. DivCG[:,:] = (vC1[:,:] + vC2[:,:]) / JC[:,:]
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        Div[ind] += DivCG[iP,jP] / CG.M[1,ind]
      end
    end
  end

  # Hyperdiffusion 
  @inbounds for iF = 1:NF
    @views JC = Global.Metric.JC[:,:,1,iF];
    @views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        DivCG[iP,jP] = Div[ind]
      end
    end
    @views mul!(D1cCG[:,:],CG.DS,DivCG[:,:])
    @views mul!(D2cCG[:,:],DivCG[:,:],CG.DST)
  
    @views @. grad1CG[:,:] = (dXdxIC[:,:,1,1,1] * D1cCG[:,:] + dXdxIC[:,:,1,2,1] * D2cCG[:,:])
    @views @. grad2CG[:,:] = (dXdxIC[:,:,1,1,2] * D1cCG[:,:] + dXdxIC[:,:,1,2,2] * D2cCG[:,:])

    @views @. D1gradCG[:,:] = dXdxIC[:,:,1,1,1] * grad1CG[:,:] + dXdxIC[:,:,1,1,2] * grad2CG[:,:]
    @views @. D2gradCG[:,:] = dXdxIC[:,:,1,2,1] * grad1CG[:,:] + dXdxIC[:,:,1,2,2] * grad2CG[:,:]

    @views mul!(vC1[:,:],CG.DW,D1gradCG[:,:])
    @views mul!(vC2[:,:],D2gradCG[:,:],CG.DWT)
    @views @. DivCG[:,:] -= HyperDDiv*(vC1[:,:] + vC2[:,:]) / JC[:,:]

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
  JF = Global.Metric.JF
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
    @views JC = Global.Metric.JC[:,:,1,iF];
    @views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]

    @views mul!(D1cCG[:,:],CG.DS,hCG[:,:,iF])
    @views mul!(D2cCG[:,:],hCG[:,:,iF],CG.DST)

    @views @. grad1CG[:,:] = dXdxIC[:,:,1,1,1] * D1cCG[:,:] + dXdxIC[:,:,1,2,1] * D2cCG[:,:]
    @views @. grad2CG[:,:] = dXdxIC[:,:,1,1,2] * D1cCG[:,:] + dXdxIC[:,:,1,2,2] * D2cCG[:,:]

    @views @. D1gradCG[:,:] = dXdxIC[:,:,1,1,1] * grad1CG[:,:] + dXdxIC[:,:,1,1,2] * grad2CG[:,:]
    @views @. D2gradCG[:,:] = dXdxIC[:,:,1,2,1] * grad1CG[:,:] + dXdxIC[:,:,1,2,2] * grad2CG[:,:]

    @views mul!(vC1[:,:],CG.DW,D1gradCG[:,:])
    @views mul!(vC2[:,:],D2gradCG[:,:],CG.DWT)
    @views @. DivCG[:,:] = HyperDDiv*(vC1[:,:] + vC2[:,:]) / JC[:,:]
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
