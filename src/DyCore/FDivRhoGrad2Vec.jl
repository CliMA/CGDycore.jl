function FDivRhoGrad2Vec(cCG,RhoCG,CG,Param)
nz=Param.Grid.nz;
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
dXdxIC = Param.cache.dXdxIC
D1cCG=reshape(
  CG.DS*reshape(cCG,OP,OP*NF*nz)
  ,OP,OP,NF,nz);
D2cCG=permute(reshape(
  CG.DS*reshape(
  permute(
  reshape(cCG,OP,OP,NF,nz)
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]);
gradCG1=RhoCG .*
  (dXdxIC[:,:,:,:,1,1].*D1cCG +
  dXdxIC[:,:,:,:,2,1].*D2cCG);
gradCG2=RhoCG .*
  (dXdxIC[:,:,:,:,1,2].*D1cCG +
  dXdxIC[:,:,:,:,2,2].*D2cCG);
divCG=(reshape(
  CG.DW*(reshape(gradCG1,OP,OP*NF*nz) .*
  reshape(dXdxIC[:,:,:,:,1,1],OP,OP*NF*nz) +
  reshape(gradCG2,OP,OP*NF*nz) .*
  reshape(dXdxIC[:,:,:,:,1,2],OP,OP*NF*nz))
  ,OP,OP,NF,nz) +
  permute(
  reshape(
  CG.DW*reshape(
  permute(
  (reshape(gradCG1,OP,OP,NF,nz) .*
  dXdxIC[:,:,:,:,2,1] +
  reshape(gradCG2,OP,OP,NF,nz) .*
  dXdxIC[:,:,:,:,2,2])
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]))./Param.JC;

return divCG
end
