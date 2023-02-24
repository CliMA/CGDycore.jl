function GradDiv!(grad1CG,grad2CG,v1CG,v2CG,Fe,Cache)
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
