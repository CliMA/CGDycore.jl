function FVort2VecDSS!(VortCG,v1CG,v2CG,CG,Global,iF)
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;

vCon = Global.Cache.CacheC1
DvCon = Global.Cache.CacheC2
@views JC = Global.Metric.JC[:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]

@inbounds for iz=1:nz
  @views @. vCon[:,:,iz] = v2CG[:,:,iz] * dXdxIC[:,:,iz,1,1] - v1CG[:,:,iz] * dXdxIC[:,:,iz,1,2]
  @views mul!(VortCG[:,:,iz],CG.DS,vCon[:,:,iz])
  @views @. vCon[:,:,iz] = v2CG[:,:,iz] * dXdxIC[:,:,iz,2,1] - v1CG[:,:,iz] * dXdxIC[:,:,iz,2,2]
  @views mul!(DvCon[:,:,iz],vCon[:,:,iz],CG.DST)
  @views @. VortCG[:,:,iz] += DvCon[:,:,iz]
end
end


function FVort2Vec!(Vort,U,CG,Global)

OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;
uPos=Global.Model.uPos
vPos=Global.Model.vPos
v1CG = Global.Cache.v1CG
v2CG = Global.Cache.v2CG
VortCG = Global.Cache.DivCG
Vort .= 0.0
@inbounds for iF=1:NF
  iG=0
  @inbounds for jP=1:OP
    @inbounds for iP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
      @inbounds for iz=1:nz
        v1CG[iP,jP,iz] = U[iz,ind,uPos]
        v2CG[iP,jP,iz] = U[iz,ind,vPos]
      end
    end
  end
  FVort2VecDSS!(VortCG,v1CG,v2CG,CG,Global,iF)

  iG=0
  @inbounds for jP=1:OP
    @inbounds for iP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
      @inbounds for iz=1:nz
        Vort[iz,ind] += VortCG[iP,jP,iz]
      end
    end
  end
end
@. Vort = Vort / CG.M
end

