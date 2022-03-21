function FRotCurl2Vec(v1CG,v2CG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
dXdxIC = Param.cache.dXdxIC
JC = Param.cache.JC

vC1=reshape(
  CG.DS*reshape(v2CG.*dXdxIC[:,:,:,:,1,1] -
  v1CG.*dXdxIC[:,:,:,:,1,2]
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz) -
  permute(
  reshape(
  CG.DS*reshape(
  permute(-v2CG.*dXdxIC[:,:,:,:,2,1] +
  v1CG.*dXdxIC[:,:,:,:,2,2]
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
RotCG[:,:,:,:,2]=(-dXdxIC[:,:,:,:,1,1].*D1cCG -
  dXdxIC[:,:,:,:,2,1].*D2cCG)./JC;
RotCG[:,:,:,:,1]=(dXdxIC[:,:,:,:,1,2].*D1cCG +
  dXdxIC[:,:,:,:,2,2].*D2cCG)./JC;
return RotCG
end
