function FGradDiv2Vec!(F,v1CG,v2CG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;

vCon = Param.CacheC1
DvCon = Param.CacheC2
vC1 = Param.CacheC3
D1vC1 = vCon
D2vC1 = DvCon
JC = Param.cache.JC

vCon .= v1CG.*Param.dXdxIC11 .+ v2CG.*Param.dXdxIC12;
mul!(reshape(vC1,OP,OP*NF*nz),CG.DS,reshape(vCon,OP,OP*nz*NF))
vCon .= v1CG.*Param.dXdxIC21 .+ v2CG.*Param.dXdxIC22;
mul!(reshape(PermutedDimsArray(DvCon,(2,1,3,4)),OP,OP*NF*nz),CG.DS,reshape(PermutedDimsArray(vCon,(2,1,3,4)),OP,OP*nz*NF))
vC1 .= vC1 .+ DvCon

mul!(reshape(D1vC1,OP,OP*NF*nz),CG.DW,reshape(vC1,OP,OP*nz*NF))
mul!(reshape(PermutedDimsArray(D2vC1,(2,1,3,4)),OP,OP*NF*nz),CG.DW,reshape(PermutedDimsArray(vC1,(2,1,3,4)),OP,OP*nz*NF))


@views F[:,:,:,:,1] .-= Param.HyperDGrad .* (Param.dXdxIC11.*D1vC1 .+ Param.dXdxIC21.*D2vC1)./JC;
@views F[:,:,:,:,2] .-= Param.HyperDGrad .* (Param.dXdxIC12.*D1vC1 .+ Param.dXdxIC22.*D2vC1)./JC;

end

function FGradDiv2Vec(v1CG,v2CG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
dXdxIC = Param.cache.dXdxIC
JC = Param.cache.JC
vC1=reshape(
  CG.DS*(reshape(v1CG.*dXdxIC[:,:,:,:,1,1] +
  v2CG.*dXdxIC[:,:,:,:,1,2]
  ,OP,OP*NF*nz))
  ,OP,OP,NF,nz) +
  permute(
  reshape(
  CG.DS*reshape(
  permute(
   v1CG.*dXdxIC[:,:,:,:,2,1] +
   v2CG.*dXdxIC[:,:,:,:,2,2]
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]);


D1cCG=reshape(
  CG.DW*reshape(vC1,OP,OP*NF*nz)
  ,OP,OP,NF,nz);
D2cCG=permute(reshape(
  CG.DW*reshape(
  permute(
  reshape(vC1,OP,OP,NF,nz)
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]);

gradCG=zeros(OP,OP,NF,nz,2);
gradCG[:,:,:,:,1]=(dXdxIC[:,:,:,:,1,1].*D1cCG +
  dXdxIC[:,:,:,:,2,1].*D2cCG)./JC;
gradCG[:,:,:,:,2]=(dXdxIC[:,:,:,:,1,2].*D1cCG +
  dXdxIC[:,:,:,:,2,2].*D2cCG)./JC;
return gradCG
end
