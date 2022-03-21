function FGradDiv2VecDSS(v1CG,v2CG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
dXdxIC = Param.cache.dXdxIC
JC = Param.cache.JC
vC1=reshape(
  CG.DS*reshape(v1CG.*dXdxIC[:,:,:,:,1,1] +
  v2CG.*dXdxIC[:,:,:,:,1,2]
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz) +
  permute(
  reshape(
  CG.DS*reshape(
  permute(v1CG.*dXdxIC[:,:,:,:,2,1] +
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
  permute(vC1
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]);

gradCG=zeros(OP,OP,NF,nz,2);
gradCG[:,:,:,:,1]=(dXdxIC[:,:,:,:,1,1].*D1cCG +
  dXdxIC[:,:,:,:,2,1].*D2cCG)./JC;
gradCG[:,:,:,:,2]=(dXdxIC[:,:,:,:,1,2].*D1cCG +
  dXdxIC[:,:,:,:,2,2].*D2cCG)./JC;


grad=zeros(CG.NumG,nz,2);
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:],OP*OP*size(i,1))
  grad[arr,:,:] = grad[arr,:,:] .+
    reshape(gradCG[:,:,i,:,:],OP*OP*size(i,1),nz,2);
end
grad[:,:,1]=grad[:,:,1]./CG.M;
grad[:,:,2]=grad[:,:,2]./CG.M;
grad1CG=reshape(grad[reshape(CG.Glob,OP*OP*NF,1),:,1]
  ,OP,OP,NF,nz);
grad2CG=reshape(grad[reshape(CG.Glob,OP*OP*NF,1),:,2]
  ,OP,OP,NF,nz);
return (grad1CG,grad2CG)
end
