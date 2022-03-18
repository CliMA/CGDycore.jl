function FVort2VecDSS(v1,v2,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
dXdxIC = Param.cache.dXdxIC
v1CG=reshape(v1[reshape(CG.Glob,OP*OP*NF,1),:]
  ,OP,OP,NF,nz);
v2CG=reshape(v2[reshape(CG.Glob,OP*OP*NF,1),:]
  ,OP,OP,NF,nz);
vC1=reshape(
  CG.DS*reshape(dXdxIC[:,:,:,:,1,1].*v2CG -
  dXdxIC[:,:,:,:,1,2].*v1CG
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

Vort=zeros(CG.NumG,nz);
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:],OP*OP*size(i,1))
  Vort[arr,:] = Vort[arr,:] .+
    reshape(vC1[:,:,i,:]
    ,OP*OP*size(i,1),nz);
end
Vort=Vort./CG.M;
return Vort
end
