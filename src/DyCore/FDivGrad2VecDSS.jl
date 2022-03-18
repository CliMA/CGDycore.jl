function FDivGrad2VecDSS(cCG,CG,Param)
nz=Param.Grid.nz;
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
D1cCG=reshape(
  CG.DS*reshape(cCG,OP,OP*NF*nz)
  ,OP,OP,NF,nz);
D2cCG=permute(reshape(
  CG.DS
  *reshape(
  permute(
  reshape(cCG,OP,OP,NF,nz)
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]);
gradCG1=Param.dXdxIC[:,:,:,:,1,1].*D1cCG +
  Param.dXdxIC[:,:,:,:,2,1].*D2cCG;
gradCG2=Param.dXdxIC[:,:,:,:,1,2].*D1cCG +
  Param.dXdxIC[:,:,:,:,2,2].*D2cCG;
vC1=(reshape(
  CG.DW*(reshape(gradCG1,OP,OP*NF*nz) .*
  reshape(Param.dXdxIC[:,:,:,:,1,1],OP,OP*NF*nz) +
  reshape(gradCG2,OP,OP*NF*nz) .*
  reshape(Param.dXdxIC[:,:,:,:,1,2],OP,OP*NF*nz))
  ,OP,OP,NF,nz) +
  permute(
  reshape(
  CG.DW*reshape(
  permute(
  (reshape(gradCG1,OP,OP,NF,nz).*Param.dXdxIC[:,:,:,:,2,1] +
   reshape(gradCG2,OP,OP,NF,nz).*Param.dXdxIC[:,:,:,:,2,2])
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]))./Param.JC;
  
  
div=zeros(CG.NumG,nz);
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:],OP*OP*size(i,1))
  div[arr,:] = div[arr,:] .+
    reshape(vC1[:,:,i,:]
    ,OP*OP*size(i,1),nz);
end
div=div./CG.M;
divCG=reshape(div[reshape(CG.Glob,OP*OP*NF,1),:]
  ,OP,OP,NF,nz);
return divCG
end
