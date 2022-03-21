function FRotCurl2Vec!(F,v1CG,v2CG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;

vCon = Param.CacheC1
DvCon = Param.CacheC2
vC1 = Param.CacheC3
D1vC1 = vCon
D2vC1 = DvCon

vCon .= v2CG.*Param.dXdxIC11 .- v1CG.*Param.dXdxIC12
mul!(reshape(vC1,OP,OP*NF*nz),CG.DS,reshape(vCon,OP,OP*nz*NF))
vCon .=  .-v2CG.*Param.dXdxIC21 .+ v1CG.*Param.dXdxIC22
mul!(reshape(PermutedDimsArray(DvCon,(2,1,3,4)),OP,OP*NF*nz),CG.DS,reshape(PermutedDimsArray(vCon,(2,1,3,4)),OP,OP*nz*NF))
vC1 .= vC1 .+ DvCon

mul!(reshape(D1vC1,OP,OP*NF*nz),CG.DW,reshape(vC1,OP,OP*nz*NF))
mul!(reshape(PermutedDimsArray(D2vC1,(2,1,3,4)),OP,OP*NF*nz),CG.DS,reshape(PermutedDimsArray(vCon,(2,1,3,4)),OP,OP*nz*NF))
@views F(:,:,:,:,1) .= (.-Param.dXdxIC11.*D1cCG .-
  Param.dXdxIC21.*D2cCG)./Param.JC;
@views F(:,:,:,:,2) .= (Param.dXdxIC12.*D1cCG +
  Param.dXdxIC22.*D2cCG)./Param.JC;

end

function FRotCurl2Vec(v1CG,v2CG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;

vC1=reshape(
  CG.DS*reshape(v2CG.*Param.dXdxIC[:,:,:,:,1,1] -
  v1CG.*Param.dXdxIC[:,:,:,:,1,2]
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz) -
  permute(
  reshape(
  CG.DS*reshape(
  permute(-v2CG.*Param.dXdxIC[:,:,:,:,2,1] +
  v1CG.*Param.dXdxIC[:,:,:,:,2,2]
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]);

D1cCG=reshape(
  CG.DW*reshape(vC1,OP,OP*NF*nz)
  ,OP,OP,NF,nz);
D2cCG=permute(reshape(
  CG.DW
  *reshape(
  permute(
  reshape(vC1,OP,OP,NF,nz)
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]);
RotCG=zeros(OP,OP,NF,nz,2);
RotCG[:,:,:,:,2]=(-Param.dXdxIC[:,:,:,:,1,1].*D1cCG -
  Param.dXdxIC[:,:,:,:,2,1].*D2cCG)./Param.JC;
RotCG[:,:,:,:,1]=(Param.dXdxIC[:,:,:,:,1,2].*D1cCG +
  Param.dXdxIC[:,:,:,:,2,2].*D2cCG)./Param.JC;
end

