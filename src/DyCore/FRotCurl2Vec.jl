function FRotCurl2Vec!(F,v1CG,v2CG,CG,Global,iF)
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;

vCon = Global.Cache.CacheC1
DvCon = Global.Cache.CacheC2
vC1 = Global.Cache.CacheC3
D1vC1 = vCon
D2vC1 = DvCon
@views JC = Global.Metric.JC[:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]

@inbounds for iz=1:nz
  @views @. vCon[:,:,iz] = v2CG[:,:,iz] * dXdxIC[:,:,iz,1,1] - v1CG[:,:,iz] * dXdxIC[:,:,iz,1,2]
  @views mul!(vC1[:,:,iz],CG.DS,vCon[:,:,iz])
  @views @. vCon[:,:,iz] = v2CG[:,:,iz] * dXdxIC[:,:,iz,2,1] - v1CG[:,:,iz] * dXdxIC[:,:,iz,2,2]
  @views mul!(DvCon[:,:,iz],vCon[:,:,iz],CG.DST)
  @views @. vC1[:,:,iz] = vC1[:,:,iz] + DvCon[:,:,iz]

  @views mul!(D1vC1[:,:,iz],CG.DW,vC1[:,:,iz])
  @views mul!(D2vC1[:,:,iz],vC1[:,:,iz],CG.DWT)
  @views @. F[:,:,iz,2] -= Global.Model.HyperDCurl*(-dXdxIC[:,:,iz,1,1] * D1vC1[:,:,iz] -
    dXdxIC[:,:,iz,2,1] * D2vC1[:,:,iz]) / JC[:,:,iz];
  @views @. F[:,:,iz,1] -= Global.Model.HyperDCurl*(dXdxIC[:,:,iz,1,2] * D1vC1[:,:,iz] +
    dXdxIC[:,:,iz,2,2] * D2vC1[:,:,iz]) / JC[:,:,iz];
end
end  

function FRotCurl2Vec1!(F,v1CG,v2CG,CG,Global)
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;

vCon = Global.Cache.CacheC1
DvCon = Global.Cache.CacheC2
vC1 = Global.Cache.CacheC3
D1vC1 = vCon
D2vC1 = DvCon
JC = Global.Metric.JC 
dXdxIC = Global.Metric.dXdxIC

@inbounds for iF=1:NF
  @inbounds for iz=1:nz
    @views @. vCon[:,:,iz,iF] = v2CG[:,:,iz,iF] * dXdxIC[:,:,iz,1,1,iF] - v1CG[:,:,iz,iF] * dXdxIC[:,:,iz,1,2,iF]
    @views mul!(vC1[:,:,iz,iF],CG.DS,vCon[:,:,iz,iF])
    @views @. vCon[:,:,iz,iF] = v2CG[:,:,iz,iF] * dXdxIC[:,:,iz,2,1,iF] - v1CG[:,:,iz,iF] * dXdxIC[:,:,iz,2,2,iF]
    @views mul!(DvCon[:,:,iz,iF],vCon[:,:,iz,iF],CG.DST)
    @views @. vC1[:,:,iz,iF] = vC1[:,:,iz,iF] + DvCon[:,:,iz,iF]

    @views mul!(D1vC1[:,:,iz,iF],CG.DW,vC1[:,:,iz,iF])
    @views mul!(D2vC1[:,:,iz,iF],vC1[:,:,iz,iF],CG.DWT)
    @views @. F[:,:,iz,iF,2] -= Global.Model.HyperDCurl*(-dXdxIC[:,:,iz,1,1,iF] * D1vC1[:,:,iz,iF] -
      dXdxIC[:,:,iz,2,1,iF] * D2vC1[:,:,iz,iF]) / JC[:,:,iz,iF];
    @views @. F[:,:,iz,iF,1] -= Global.Model.HyperDCurl*(dXdxIC[:,:,iz,1,2,iF] * D1vC1[:,:,iz,iF] +
      dXdxIC[:,:,iz,2,2,iF] * D2vC1[:,:,iz,iF]) / JC[:,:,iz,iF];
  end
end  
end

