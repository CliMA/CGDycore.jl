function FDivRhoGrad2Vec!(F,cCG,RhoCG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = Global.ThreadCache
nz=Global.Grid.nz;
OP=CG.OrdPoly+1;

D1cCG = TCacheC1[Threads.threadid()]
D2cCG = TCacheC2[Threads.threadid()]
grad1CG = TCacheC3[Threads.threadid()]
grad2CG = TCacheC4[Threads.threadid()]
D1gradCG = TCacheC1[Threads.threadid()]
D2gradCG = TCacheC2[Threads.threadid()]
vC1 = TCacheC3[Threads.threadid()]
vC2 = TCacheC4[Threads.threadid()]

@views JC = Global.Metric.JC[:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]

@inbounds for iz=1:nz
  @views mul!(D1cCG[:,:],CG.DS,cCG[:,:,iz])
  @views mul!(D2cCG[:,:],cCG[:,:,iz],CG.DST)

  @views @. grad1CG[:,:] = RhoCG[:,:,iz] * (dXdxIC[:,:,iz,1,1] * D1cCG[:,:] + dXdxIC[:,:,iz,2,1] * D2cCG[:,:])
  @views @. grad2CG[:,:] = RhoCG[:,:,iz] * (dXdxIC[:,:,iz,1,2] * D1cCG[:,:] + dXdxIC[:,:,iz,2,2] * D2cCG[:,:])

  @views @. D1gradCG[:,:] = dXdxIC[:,:,iz,1,1] * grad1CG[:,:] + dXdxIC[:,:,iz,1,2] * grad2CG[:,:]
  @views @. D2gradCG[:,:] = dXdxIC[:,:,iz,2,1] * grad1CG[:,:] + dXdxIC[:,:,iz,2,2] * grad2CG[:,:]

  @views mul!(vC1[:,:],CG.DW,D1gradCG[:,:])
  @views mul!(vC2[:,:],D2gradCG[:,:],CG.DWT)
  @views @. F[:,:,iz] -= Global.Model.HyperDDiv*(vC1[:,:] + vC2[:,:]) / JC[:,:,iz]

end
end

function FDivGrad2VecDSS!(divCG,cCG,RhoCG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4, TCacheC5 = Global.ThreadCache
nz=Global.Grid.nz;
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
D1cCG = TCacheC1[Threads.threadid()]
D2cCG = TCacheC2[Threads.threadid()]
grad1CG = TCacheC3[Threads.threadid()]
grad2CG = TCacheC4[Threads.threadid()]
ccCG = TCacheC5[Threads.threadid()]
D1gradCG = TCacheC1[Threads.threadid()]
D2gradCG = TCacheC2[Threads.threadid()]
vC1 = TCacheC3[Threads.threadid()]
vC2 = TCacheC4[Threads.threadid()]
@views JC = Global.Metric.JC[:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]

@inbounds for iz=1:nz
  @views @. ccCG[:,:] = cCG[:,:,iz] / RhoCG[:,:,iz]
  @views mul!(D1cCG[:,:],CG.DS,ccCG[:,:])
  @views mul!(D2cCG[:,:],ccCG[:,:],CG.DST)

  @views @. grad1CG[:,:] = dXdxIC[:,:,iz,1,1] * D1cCG[:,:] + dXdxIC[:,:,iz,2,1] * D2cCG[:,:]
  @views @. grad2CG[:,:] = dXdxIC[:,:,iz,1,2] * D1cCG[:,:] + dXdxIC[:,:,iz,2,2] * D2cCG[:,:]

  @views @. D1gradCG[:,:] = dXdxIC[:,:,iz,1,1] * grad1CG[:,:] + dXdxIC[:,:,iz,1,2] * grad2CG[:,:]
  @views @. D2gradCG[:,:] = dXdxIC[:,:,iz,2,1] * grad1CG[:,:] + dXdxIC[:,:,iz,2,2] * grad2CG[:,:]

  @views mul!(vC1[:,:],CG.DW,D1gradCG[:,:])
  @views mul!(vC2[:,:],D2gradCG[:,:],CG.DWT)
  @views @. divCG[:,:,iz] = (vC1[:,:] + vC2[:,:]) / JC[:,:,iz]
end
end

function FDivGrad2VecDSS!(divCG,cCG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4, TCacheC5 = Global.ThreadCache
nz=Global.Grid.nz;
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
D1cCG = TCacheC1[Threads.threadid()]
D2cCG = TCacheC2[Threads.threadid()]
grad1CG = TCacheC3[Threads.threadid()]
grad2CG = TCacheC4[Threads.threadid()]
ccCG = TCacheC5[Threads.threadid()]
D1gradCG = TCacheC1[Threads.threadid()]
D2gradCG = TCacheC2[Threads.threadid()]
vC1 = TCacheC3[Threads.threadid()]
vC2 = TCacheC4[Threads.threadid()]
@views JC = Global.Metric.JC[:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]

@inbounds for iz=1:nz
  @views @. ccCG[:,:] = cCG[:,:,iz] 
  @views mul!(D1cCG[:,:],CG.DS,ccCG[:,:])
  @views mul!(D2cCG[:,:],ccCG[:,:],CG.DST)

  @views @. grad1CG[:,:] = dXdxIC[:,:,iz,1,1] * D1cCG[:,:] + dXdxIC[:,:,iz,2,1] * D2cCG[:,:]
  @views @. grad2CG[:,:] = dXdxIC[:,:,iz,1,2] * D1cCG[:,:] + dXdxIC[:,:,iz,2,2] * D2cCG[:,:]

  @views @. D1gradCG[:,:] = dXdxIC[:,:,iz,1,1] * grad1CG[:,:] + dXdxIC[:,:,iz,1,2] * grad2CG[:,:]
  @views @. D2gradCG[:,:] = dXdxIC[:,:,iz,2,1] * grad1CG[:,:] + dXdxIC[:,:,iz,2,2] * grad2CG[:,:]

  @views mul!(vC1[:,:],CG.DW,D1gradCG[:,:])
  @views mul!(vC2[:,:],D2gradCG[:,:],CG.DWT)
  @show size(divCG)
  @show nz
  @show iz
  @views @. divCG[:,:,iz] = (vC1[:,:] + vC2[:,:]) / JC[:,:,iz]
end
end


function FDivGrad2VecDSS!(divCG,cCG,RhoCG,dXdxIC,JC,CG,ThreadCache)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4, TCacheC5 = ThreadCache

D1cCG = TCacheC1[Threads.threadid()]
D2cCG = TCacheC2[Threads.threadid()]
grad1CG = TCacheC3[Threads.threadid()]
grad2CG = TCacheC4[Threads.threadid()]
ccCG = TCacheC5[Threads.threadid()]
D1gradCG = TCacheC1[Threads.threadid()]
D2gradCG = TCacheC2[Threads.threadid()]
vC1 = TCacheC3[Threads.threadid()]
vC2 = TCacheC4[Threads.threadid()]

  @. @fastmath ccCG = cCG / RhoCG
  mul!(D1cCG,CG.DS,ccCG)
  mul!(D2cCG,ccCG,CG.DST)

  @views @fastmath @. grad1CG = dXdxIC[:,:,1,1] * D1cCG + dXdxIC[:,:,2,1] * D2cCG[:,:]
  @views @fastmath @. grad2CG = dXdxIC[:,:,1,2] * D1cCG + dXdxIC[:,:,2,2] * D2cCG[:,:]

  @views @fastmath @. D1gradCG = dXdxIC[:,:,1,1] * grad1CG + dXdxIC[:,:,1,2] * grad2CG[:,:]
  @views @fastmath @. D2gradCG = dXdxIC[:,:,2,1] * grad1CG + dXdxIC[:,:,2,2] * grad2CG[:,:]

  mul!(vC1,CG.DW,D1gradCG)
  mul!(vC2,D2gradCG,CG.DWT)
  @fastmath @. divCG = (vC1 + vC2) / JC
end


