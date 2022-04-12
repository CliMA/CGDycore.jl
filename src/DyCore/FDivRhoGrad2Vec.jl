function FDivRhoGrad2Vec!(F,cCG,RhoCG,CG,Global,iF)
nz=Global.Grid.nz;
OP=CG.OrdPoly+1;

D1cCG = Global.Cache.CacheC1
D2cCG = Global.Cache.CacheC2
grad1CG = Global.Cache.CacheC3
grad2CG = Global.Cache.CacheC4
D1gradCG = Global.Cache.CacheC1
D2gradCG = Global.Cache.CacheC2
vC1 = Global.Cache.CacheC3
vC2 = Global.Cache.CacheC4
@views JC = Global.Metric.JC[:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]


@inbounds for iz=1:nz
  @views mul!(D1cCG[:,:,iz],CG.DS,cCG[:,:,iz])
  @views mul!(D2cCG[:,:,iz],cCG[:,:,iz],CG.DST)

  @views @. grad1CG[:,:,iz] = RhoCG[:,:,iz] * (dXdxIC[:,:,iz,1,1] * D1cCG[:,:,iz] + dXdxIC[:,:,iz,2,1] * D2cCG[:,:,iz])
  @views @. grad2CG[:,:,iz] = RhoCG[:,:,iz] * (dXdxIC[:,:,iz,1,2] * D1cCG[:,:,iz] + dXdxIC[:,:,iz,2,2] * D2cCG[:,:,iz])

  @views @. D1gradCG[:,:,iz] = dXdxIC[:,:,iz,1,1] * grad1CG[:,:,iz] + dXdxIC[:,:,iz,1,2] * grad2CG[:,:,iz]
  @views @. D2gradCG[:,:,iz] = dXdxIC[:,:,iz,2,1] * grad1CG[:,:,iz] + dXdxIC[:,:,iz,2,2] * grad2CG[:,:,iz]

  @views mul!(vC1[:,:,iz],CG.DW,D1gradCG[:,:,iz])
  @views mul!(vC2[:,:,iz],D2gradCG[:,:,iz],CG.DWT)
  @views @. F[:,:,iz] -= Global.Model.HyperDDiv*(vC1[:,:,iz] + vC2[:,:,iz]) / JC[:,:,iz]

end

end


