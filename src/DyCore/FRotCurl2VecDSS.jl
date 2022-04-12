function FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;

vCon = Global.Cache.CacheC1
DvCon = Global.Cache.CacheC2
vC1 = Global.Cache.CacheC3
D1vC1 = vCon
D2vC1 = DvCon
Rot1TempCG = Global.Cache.CacheC3
Rot2TempCG = Global.Cache.CacheC4
Rot1 = Global.Cache.Cache1
Rot2 = Global.Cache.Cache2
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
  @views @. Rot2CG[:,:,iz] = (-dXdxIC[:,:,iz,1,1] * D1vC1[:,:,iz] - 
    dXdxIC[:,:,iz,2,1] * D2vC1[:,:,iz]) / JC[:,:,iz];
  @views @. Rot1CG[:,:,iz] = (dXdxIC[:,:,iz,1,2] * D1vC1[:,:,iz] + 
    dXdxIC[:,:,iz,2,2] * D2vC1[:,:,iz]) / JC[:,:,iz];
end

end


function FRotCurl2VecDSS1!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global)
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;

vCon = Global.Cache.CacheC1
DvCon = Global.Cache.CacheC2
vC1 = Global.Cache.CacheC3
D1vC1 = vCon
D2vC1 = DvCon
Rot1TempCG = Global.Cache.CacheC3
Rot2TempCG = Global.Cache.CacheC4
Rot1 = Global.Cache.Cache1
Rot2 = Global.Cache.Cache2
JC = Global.Metric.JC;
dXdxIC = Global.Metric.dXdxIC

@. Rot1 = 0
@. Rot2 = 0
@inbounds for iF=1:NF
  @inbounds for iz=1:nz  
    @views @. vCon[:,:,iz,iF] = v2CG[:,:,iz,iF] * dXdxIC[:,:,iz,1,1,iF] - v1CG[:,:,iz,iF] * dXdxIC[:,:,iz,1,2,iF]
    @views mul!(vC1[:,:,iz,iF],CG.DS,vCon[:,:,iz,iF])
    @views @. vCon[:,:,iz,iF] = v2CG[:,:,iz,iF] * dXdxIC[:,:,iz,2,1,iF] - v1CG[:,:,iz,iF] * dXdxIC[:,:,iz,2,2,iF]
    @views mul!(DvCon[:,:,iz,iF],vCon[:,:,iz,iF],CG.DST)
    @views @. vC1[:,:,iz,iF] = vC1[:,:,iz,iF] + DvCon[:,:,iz,iF]

    @views mul!(D1vC1[:,:,iz,iF],CG.DW,vC1[:,:,iz,iF])
    @views mul!(D2vC1[:,:,iz,iF],vC1[:,:,iz,iF],CG.DWT)
    @views @. Rot2TempCG[:,:,iz,iF] = (-dXdxIC[:,:,iz,1,1,iF] * D1vC1[:,:,iz,iF] - 
      dXdxIC[:,:,iz,2,1,iF] * D2vC1[:,:,iz,iF]) / JC[:,:,iz,iF];
    @views @. Rot1TempCG[:,:,iz,iF] = (dXdxIC[:,:,iz,1,2,iF] * D1vC1[:,:,iz,iF] + 
      dXdxIC[:,:,iz,2,2,iF] * D2vC1[:,:,iz,iF]) / JC[:,:,iz,iF];
  end

  iG=0
  @inbounds for iP=1:OP
    @inbounds for jP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
      @inbounds for iz=1:nz
        Rot1[iz,ind] += Rot1TempCG[iP,jP,iz,iF]
        Rot2[iz,ind] += Rot2TempCG[iP,jP,iz,iF]
      end
    end
  end
end

@. Rot1 = Rot1 / CG.M;
@. Rot2 = Rot2 / CG.M;
@inbounds for iF=1:NF
  iG=0
  @inbounds for iP=1:OP
    @inbounds for jP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
      @inbounds for iz=1:nz
        Rot1CG[iP,jP,iz,iF] = Rot1[iG,iz]
        Rot2CG[iP,jP,iz,iF] = Rot2[iG,iz]
      end
    end
  end
end
end

