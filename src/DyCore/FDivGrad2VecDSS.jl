function FDivGrad2VecDSS!(divCG,cCG,RhoCG,CG,Global,iF)
nz=Global.Grid.nz;
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;

D1cCG = Global.Cache.CacheC1
D2cCG = Global.Cache.CacheC2
grad1CG = Global.Cache.CacheC3
grad2CG = Global.Cache.CacheC4
ccCG = Global.Cache.CacheC5
D1gradCG = Global.Cache.CacheC1
D2gradCG = Global.Cache.CacheC2
vC1 = Global.Cache.CacheC3
vC2 = Global.Cache.CacheC4
@views JC = Global.Metric.JC[:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]

@inbounds for iz=1:nz  
  @views @. ccCG[:,:,iz] = cCG[:,:,iz] / RhoCG[:,:,iz]   
  @views mul!(D1cCG[:,:,iz],CG.DS,ccCG[:,:,iz])
  @views mul!(D2cCG[:,:,iz],ccCG[:,:,iz],CG.DST)

  @views @. grad1CG[:,:,iz] = dXdxIC[:,:,iz,1,1] * D1cCG[:,:,iz] + dXdxIC[:,:,iz,2,1] * D2cCG[:,:,iz]
  @views @. grad2CG[:,:,iz] = dXdxIC[:,:,iz,1,2] * D1cCG[:,:,iz] + dXdxIC[:,:,iz,2,2] * D2cCG[:,:,iz]

  @views @. D1gradCG[:,:,iz] = dXdxIC[:,:,iz,1,1] * grad1CG[:,:,iz] + dXdxIC[:,:,iz,1,2] * grad2CG[:,:,iz]
  @views @. D2gradCG[:,:,iz] = dXdxIC[:,:,iz,2,1] * grad1CG[:,:,iz] + dXdxIC[:,:,iz,2,2] * grad2CG[:,:,iz]

  @views mul!(vC1[:,:,iz],CG.DW,D1gradCG[:,:,iz])
  @views mul!(vC2[:,:,iz],D2gradCG[:,:,iz],CG.DWT)
  @views @. divCG[:,:,iz] = (vC1[:,:,iz] + vC2[:,:,iz]) / JC[:,:,iz]
end
end

function FDivGrad2VecDSS1!(divCG,cCG,RhoCG,CG,Global)
nz=Global.Grid.nz;
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;

D1cCG = Global.Cache.CacheC1
D2cCG = Global.Cache.CacheC2
grad1CG = Global.Cache.CacheC3
grad2CG = Global.Cache.CacheC4
ccCG = Global.Cache.CacheC5
D1gradCG = Global.Cache.CacheC1
D2gradCG = Global.Cache.CacheC2
vC1 = Global.Cache.CacheC3
vC2 = Global.Cache.CacheC4
div = Global.Cache.Cache1
JC = Global.Metric.JC;
dXdxIC = Global.Metric.dXdxIC;

@. div = 0
@inbounds for iF=1:NF
  @inbounds for iz=1:nz  
    @views @. ccCG[:,:,iz,iF] = cCG[:,:,iz,iF] / RhoCG[:,:,iz,iF]   
    @views mul!(D1cCG[:,:,iz,iF],CG.DS,ccCG[:,:,iz,iF])
    @views mul!(D2cCG[:,:,iz,iF],ccCG[:,:,iz,iF],CG.DST)

    @views @. grad1CG[:,:,iz,iF] = dXdxIC[:,:,iz,1,1,iF] * D1cCG[:,:,iz,iF] + dXdxIC[:,:,iz,2,1,iF] * D2cCG[:,:,iz,iF]
    @views @. grad2CG[:,:,iz,iF] = dXdxIC[:,:,iz,1,2,iF] * D1cCG[:,:,iz,iF] + dXdxIC[:,:,iz,2,2,iF] * D2cCG[:,:,iz,iF]

    @views @. D1gradCG[:,:,iz,iF] = dXdxIC[:,:,iz,1,1,iF] * grad1CG[:,:,iz,iF] + dXdxIC[:,:,iz,1,2,iF] * grad2CG[:,:,iz,iF]
    @views @. D2gradCG[:,:,iz,iF] = dXdxIC[:,:,iz,2,1,iF] * grad1CG[:,:,iz,iF] + dXdxIC[:,:,iz,2,2,iF] * grad2CG[:,:,iz,iF]

    @views mul!(vC1[:,:,iz,iF],CG.DW,D1gradCG[:,:,iz,iF])
    @views mul!(vC2[:,:,iz,iF],D2gradCG[:,:,iz,iF],CG.DWT)
    @views @. vC1[:,:,iz,iF] = (vC1[:,:,iz,iF] + vC2[:,:,iz,iF]) / JC[:,:,iz,iF]

  end
  iG=0
  @inbounds for iP=1:OP
    @inbounds for jP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
      @inbounds for iz=1:nz
        div[iz,ind] += vC1[iP,jP,iz,iF]
      end
    end
  end
end

@. div = div / CG.M
@inbounds for iF=1:NF
  iG=0
  @inbounds for iP=1:OP
    @inbounds for jP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
      @inbounds for iz=1:nz
        divCG[iP,jP,iz,iF] = div[iz,iG]
      end
    end
  end
end
end

