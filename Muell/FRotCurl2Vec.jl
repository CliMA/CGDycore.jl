function FRotCurl2Vec!(F,v1CG,v2CG,CG,Global,iF)
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
  @views @. vCon = v2CG[:,:,iz] * dXdxIC[:,:,iz,1,1] - v1CG[:,:,iz] * dXdxIC[:,:,iz,1,2]
  mul!(vC1,CG.DS,vCon)
  @views @. vCon = v2CG[:,:,iz] * dXdxIC[:,:,iz,2,1] - v1CG[:,:,iz] * dXdxIC[:,:,iz,2,2]
  mul!(DvCon,vCon,CG.DST)
  @views @. vC1 = vC1 + DvCon

  @views mul!(D1vC1,CG.DW,vC1)
  @views mul!(D2vC1,vC1,CG.DWT)
  @views @. F[:,:,iz,2] -= Global.Model.HyperDCurl*(-dXdxIC[:,:,iz,1,1] * D1vC1 -
    dXdxIC[:,:,iz,2,1] * D2vC1) / JC[:,:,iz];
  @views @. F[:,:,iz,1] -= Global.Model.HyperDCurl*( dXdxIC[:,:,iz,1,2] * D1vC1 +
    dXdxIC[:,:,iz,2,2] * D2vC1) / JC[:,:,iz];
end
end  

function FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
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
  @views @. vCon = v2CG[:,:,iz] * dXdxIC[:,:,iz,1,1] - v1CG[:,:,iz] * dXdxIC[:,:,iz,1,2]
  @views mul!(vC1,CG.DS,vCon)
  @views @. vCon = v2CG[:,:,iz] * dXdxIC[:,:,iz,2,1] - v1CG[:,:,iz] * dXdxIC[:,:,iz,2,2]
  @views mul!(DvCon,vCon,CG.DST)
  @views @. vC1 = vC1 + DvCon

  @views mul!(D1vC1,CG.DW,vC1)
  @views mul!(D2vC1,vC1,CG.DWT)
  @views @. Rot2CG[:,:,iz] = (-dXdxIC[:,:,iz,1,1] * D1vC1 -
    dXdxIC[:,:,iz,2,1] * D2vC1) / JC[:,:,iz];
  @views @. Rot1CG[:,:,iz] = (dXdxIC[:,:,iz,1,2] * D1vC1 +
    dXdxIC[:,:,iz,2,2] * D2vC1) / JC[:,:,iz];
end

end


function FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,dXdxIC,JC,CG,ThreadCache)
@unpack TCacheC1, TCacheC2, TCacheC3 = ThreadCache

vCon = TCacheC1[Threads.threadid()]
DvCon = TCacheC2[Threads.threadid()]
vC1 = TCacheC3[Threads.threadid()]
D1vC1 = TCacheC1[Threads.threadid()]
D2vC1 = TCacheC2[Threads.threadid()]

  @views @. vCon = v2CG * dXdxIC[:,:,1,1] - v1CG * dXdxIC[:,:,1,2]
  @views mul!(vC1,CG.DS,vCon)
  @views @. vCon = v2CG * dXdxIC[:,:,2,1] - v1CG * dXdxIC[:,:,2,2]
  @views mul!(DvCon,vCon,CG.DST)
  @views @. vC1 = vC1 + DvCon

  @views mul!(D1vC1,CG.DW,vC1)
  @views mul!(D2vC1,vC1,CG.DWT)
  @views @. Rot2CG = (-dXdxIC[:,:,1,1] * D1vC1 -
    dXdxIC[:,:,2,1] * D2vC1) / JC
  @views @. Rot1CG = (dXdxIC[:,:,1,2] * D1vC1 +
    dXdxIC[:,:,2,2] * D2vC1) / JC

end


