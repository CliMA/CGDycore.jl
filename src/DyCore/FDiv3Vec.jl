function FDiv3Vec!(F,cCG,v1CG,v2CG,v3CG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
dXdxIF = Param.cache.dXdxIF;
dXdxIC = Param.cache.dXdxIC;
# Contravariant components

vCon = Param.CacheC1
DvCon = Param.CacheC2

vCon .= (v1CG.*Param.dXdxIC11 .+ v2CG.*Param.dXdxIC12) .* cCG;
mul!(reshape(DvCon,OP,OP*NF*nz),CG.DS,reshape(vCon,OP,OP*nz*NF))
F .= F .- DvCon


vCon .= (v1CG.*Param.dXdxIC21 .+ v2CG.*Param.dXdxIC22) .* cCG;
mul!(reshape(PermutedDimsArray(DvCon,(2,1,3,4)),OP,OP*NF*nz),CG.DS,reshape(PermutedDimsArray(vCon,(2,1,3,4)),OP,OP*nz*NF))
F .= F .- DvCon

@views vCon[:,:,:,1:end-1] .= 0.5.*((v1CG[:,:,:,1:end-1] .+ v1CG[:,:,:,2:end]).*
        Param.dXdxIF[:,:,:,2:end-1,3,1] .+
      (v2CG[:,:,:,1:end-1] .+ v2CG[:,:,:,2:end]).*
        Param.dXdxIF[:,:,:,2:end-1,3,2]) .+
       v3CG[:,:,:,2:end-1].*
        Param.dXdxIF[:,:,:,2:end-1,3,3];
@views vCon[:,:,:,1:end-1] .= 0.5.*(cCG[:,:,:,1:end-1] .+ cCG[:,:,:,2:end]).*vCon[:,:,:,1:end-1];
if nz>1
  @views DvCon[:,:,:,1] .= 0.5.*vCon[:,:,:,1];
  @views DvCon[:,:,:,2:end-1] .= 0.5.*(vCon[:,:,:,2:end-1] .- vCon[:,:,:,1:end-2]);
  @views DvCon[:,:,:,end] .= -0.5.*vCon[:,:,:,end-1]; # 0.5 Metric
else
  @views DvCon[:,:,:,1] .= 0
end
F .= F .- DvCon

end

function FDiv3Vec(cCG,v1CG,v2CG,v3CG,CG,Param)
# noch uberarbeiten
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
# Contravariant components
v1Con=v1CG.*
    Param.dXdxIC[:,:,:,:,1,1] +
   v2CG.*
    Param.dXdxIC[:,:,:,:,1,2];
v2Con=v1CG.*
    Param.dXdxIC[:,:,:,:,2,1] +
   v2CG.*
    Param.dXdxIC[:,:,:,:,2,2];


v3Con=0.5*((v1CG[:,:,:,1:end-1]+v1CG[:,:,:,2:end]).*
        Param.dXdxIF[:,:,:,2:nz,3,1]+
      (v2CG[:,:,:,1:end-1]+v2CG[:,:,:,2:end]).*
        Param.dXdxIF[:,:,:,2:nz,3,2]) +
       v3CG[:,:,:,2:nz].*
        Param.dXdxIF[:,:,:,2:nz,3,3];

Dv1Con=reshape(
  CG.DS*reshape(v1Con.*cCG,OP,OP*NF*nz),
   OP,OP,NF,nz);
Dv2Con=permute(
  reshape(
  CG.DS*reshape(
  permute(v2Con.*cCG,
  [2,1,3,4]),
  OP,OP*NF*nz),
  OP,OP,NF,nz),
  [2,1,3,4]);
v3Con=0.5*(cCG[:,:,:,1:nz-1]+cCG[:,:,:,2:nz]).*v3Con;
Dv3Con=zeros(OP,OP,NF,nz);
if nz>1
  Dv3Con[:,:,:,1]=0.5*v3Con[:,:,:,1];
  Dv3Con[:,:,:,2:nz-1]=0.5*(v3Con[:,:,:,2:nz-1]-v3Con[:,:,:,1:nz-2]);
  Dv3Con[:,:,:,nz]=-0.5*v3Con[:,:,:,nz-1]; # 0.5 Metric
end
divCG=Dv1Con+Dv2Con+Dv3Con;

return divCG
end
