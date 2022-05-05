function FGradDiv2Vec!(F,v1CG,v2CG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3 = Global.ThreadCache
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;

vCon = TCacheC1[Threads.threadid()]
DvCon = TCacheC2[Threads.threadid()]
vC1 = TCacheC3[Threads.threadid()]
D1vC1 = TCacheC1[Threads.threadid()]
D2vC1 = TCacheC2[Threads.threadid()]
@views JC = Global.Metric.JC[:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]

@inbounds for iz=1:nz
  @views @. vCon = v1CG[:,:,iz] * dXdxIC[:,:,iz,1,1] + v2CG[:,:,iz] * dXdxIC[:,:,iz,1,2]
  mul!(vC1,CG.DS,vCon)
  @views @. vCon = v1CG[:,:,iz] * dXdxIC[:,:,iz,2,1] + v2CG[:,:,iz] * dXdxIC[:,:,iz,2,2]
  mul!(DvCon,vCon,CG.DST)
  @. vC1 = vC1 + DvCon

  mul!(D1vC1,CG.DW,vC1)
  mul!(D2vC1,vC1,CG.DWT)

  @views @. F[:,:,iz,1] -= Global.Model.HyperDGrad*(dXdxIC[:,:,iz,1,1] * D1vC1 +
    dXdxIC[:,:,iz,2,1] * D2vC1) / JC[:,:,iz]
  @views @. F[:,:,iz,2] -= Global.Model.HyperDGrad*(dXdxIC[:,:,iz,1,2] * D1vC1 +
    dXdxIC[:,:,iz,2,2] * D2vC1) / JC[:,:,iz]
end
end

function FGradDiv2VecDSS!(grad1CG,grad2CG,v1CG,v2CG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3 = Global.ThreadCache
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;

vCon = TCacheC1[Threads.threadid()]
DvCon = TCacheC2[Threads.threadid()]
vC1 = TCacheC3[Threads.threadid()]
D1vC1 = TCacheC1[Threads.threadid()]
D2vC1 = TCacheC2[Threads.threadid()]
@views JC = Global.Metric.JC[:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]

@inbounds for iz=1:nz
  @views @. vCon = v1CG[:,:,iz] * dXdxIC[:,:,iz,1,1] + v2CG[:,:,iz] * dXdxIC[:,:,iz,1,2]
  mul!(vC1,CG.DS,vCon)
  @views @. vCon = v1CG[:,:,iz] * dXdxIC[:,:,iz,2,1] + v2CG[:,:,iz] * dXdxIC[:,:,iz,2,2]
  mul!(DvCon,vCon,CG.DST)
  @. vC1 = vC1 + DvCon

  mul!(D1vC1,CG.DW,vC1)
  mul!(D2vC1,vC1,CG.DWT)

  @views @. grad1CG[:,:,iz] = (dXdxIC[:,:,iz,1,1] * D1vC1 +
    dXdxIC[:,:,iz,2,1] * D2vC1) / JC[:,:,iz]
  @views @. grad2CG[:,:,iz] = (dXdxIC[:,:,iz,1,2] * D1vC1 +
    dXdxIC[:,:,iz,2,2] * D2vC1) / JC[:,:,iz]
  end
end

function FGradDiv2VecDSS!(grad1CG,grad2CG,v1CG,v2CG,dXdxIC,JC,CG,ThreadCache,::Val{NPOLY}) where {NPOLY}
@unpack TCacheC1, TCacheC2, TCacheC3 = ThreadCache

  #vCon = TCacheC1[Threads.threadid()]
  DvCon = TCacheC2[Threads.threadid()]
  vC1 = TCacheC3[Threads.threadid()]
  D1vC1 = TCacheC1[Threads.threadid()]
  D2vC1 = TCacheC2[Threads.threadid()]
 
  @views vCon = SMatrix{NPOLY,NPOLY}(v1CG .* dXdxIC[:,:,1,1] .+ v2CG .* dXdxIC[:,:,1,2])
  mul!(vC1,CG.DS,vCon)
  @views vCon = SMatrix{NPOLY,NPOLY}(v1CG .* dXdxIC[:,:,2,1] .+ v2CG .* dXdxIC[:,:,2,2])
  mul!(DvCon,vCon,CG.DST)
  @. vC1 = vC1 + DvCon

  mul!(D1vC1,CG.DW,vC1)
  mul!(D2vC1,vC1,CG.DWT)

  @views @. @fastmath grad1CG = (dXdxIC[:,:,1,1] * D1vC1 +
    dXdxIC[:,:,2,1] * D2vC1) / JC
  @views @. @fastmath grad2CG = (dXdxIC[:,:,1,2] * D1vC1 +
    dXdxIC[:,:,2,2] * D2vC1) / JC
end

