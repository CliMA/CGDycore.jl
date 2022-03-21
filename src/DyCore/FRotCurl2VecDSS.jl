function FRotCurl2VecDSS(v1CG,v2CG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
dXdxIC = Param.cache.dXdxIC
JC = Param.cache.JC

vC1=reshape(
  CG.DS*reshape(dXdxIC[:,:,:,:,1,1].*v2CG -
  dXdxIC[:,:,:,:,1,2].*v1CG
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz) -
  permute(
  reshape(
  CG.DS*reshape(
  permute( -v2CG.*dXdxIC[:,:,:,:,2,1] +
  v1CG.*dXdxIC[:,:,:,:,2,2]
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]);

D1cCG=reshape(
  CG.DW*reshape(vC1,OP,OP*NF*nz)
  ,OP,OP,NF,nz);
D2cCG=permute(reshape(
  CG.DW *
  reshape(
  permute(vC1
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]);
RotCG=zeros(OP,OP,NF,nz,2);
RotCG[:,:,:,:,2]=(-dXdxIC[:,:,:,:,1,1].*D1cCG -
  dXdxIC[:,:,:,:,2,1].*D2cCG)./JC;
RotCG[:,:,:,:,1]=(dXdxIC[:,:,:,:,1,2].*D1cCG
  +dXdxIC[:,:,:,:,2,2].*D2cCG)./JC;
Rot=zeros(CG.NumG,nz,2);

for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:],OP*OP*size(i,1))
  Rot[arr,:,:] = Rot[arr,:,:] .+
    reshape(RotCG[:,:,i,:,:]
    ,OP*OP*size(i,1),nz,2);
end

Rot[:,:,1]=Rot[:,:,1]./CG.M;
Rot[:,:,2]=Rot[:,:,2]./CG.M;
Rot1CG=reshape(Rot[reshape(CG.Glob,OP*OP*NF,1),:,1]
  ,OP,OP,NF,nz);
Rot2CG=reshape(Rot[reshape(CG.Glob,OP*OP*NF,1),:,2]
  ,OP,OP,NF,nz);
return (Rot1CG,Rot2CG)
end
