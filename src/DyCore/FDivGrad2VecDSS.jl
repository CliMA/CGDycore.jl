function FDivGrad2VecDSS!(divCG,cCG,CG,Param)
nz=Param.Grid.nz;
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;

D1cCG = Param.CacheC1
D2cCG = Param.CacheC2
grad1CG = Param.CacheC3
grad2CG = Param.CacheC4
D1gradCG = Param.CacheC1
D2gradCG = Param.CacheC2
vC1 = Param.CacheC3
vC2 = Param.CacheC4
div = Param.Cache1
mul!(reshape(D1cCG,OP,OP*NF*nz),CG.DS,reshape(cCG,OP,OP*nz*NF))
mul!(reshape(PermutedDimsArray(D2cCG,(2,1,3,4)),OP,OP*NF*nz),CG.DS,reshape(PermutedDimsArray(cCG,(2,1,3,4)),OP,OP*nz*NF))

grad1CG .= Param.dXdxIC11.*D1cCG .+ Param.dXdxIC21.*D2cCG
grad2CG .= Param.dXdxIC12.*D1cCG .+ Param.dXdxIC22.*D2cCG

D1gradCG .= Param.dXdxIC11.*grad1CG .+ Param.dXdxIC12.*grad2CG
D2gradCG .= Param.dXdxIC21.*grad1CG .+ Param.dXdxIC22.*grad2CG

mul!(reshape(vC1,OP,OP*NF*nz),CG.DW,reshape(D1gradCG,OP,OP*nz*NF))
mul!(reshape(PermutedDimsArray(vC2,(2,1,3,4)),OP,OP*NF*nz),CG.DW,reshape(PermutedDimsArray(D2gradCG,(2,1,3,4)),OP,OP*nz*NF))
vC1 .= (vC1 .+ vC2) ./ Param.JC

div .= 0
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:],OP*OP*size(i,1))
  div[arr,:] .= div[arr,:] .+
    reshape(vC1[:,:,i,:]
    ,OP*OP*size(i,1),nz);
end
div=div./CG.M;
divCG .= reshape(div[reshape(CG.Glob,OP*OP*NF,1),:]
  ,OP,OP,NF,nz);
end

function FDivGrad2VecDSS(cCG,CG,Param)
nz=Param.Grid.nz;
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
dXdxIC = Param.cache.dXdxIC
JC = Param.cache.JC
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

gradCG1=dXdxIC[:,:,:,:,1,1].*D1cCG +
  dXdxIC[:,:,:,:,2,1].*D2cCG;
gradCG2=dXdxIC[:,:,:,:,1,2].*D1cCG +
  dXdxIC[:,:,:,:,2,2].*D2cCG;

vC1=(reshape(
  CG.DW*(reshape(gradCG1,OP,OP*NF*nz) .*
  reshape(dXdxIC[:,:,:,:,1,1],OP,OP*NF*nz) +
  reshape(gradCG2,OP,OP*NF*nz) .*
  reshape(dXdxIC[:,:,:,:,1,2],OP,OP*NF*nz))
  ,OP,OP,NF,nz) +
  permute(
  reshape(
  CG.DW*reshape(
  permute(
  (reshape(gradCG1,OP,OP,NF,nz).*dXdxIC[:,:,:,:,2,1] +
   reshape(gradCG2,OP,OP,NF,nz).*dXdxIC[:,:,:,:,2,2])
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]))./JC;


div=zeros(CG.NumG,nz);
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:],OP*OP*size(i,1))
  div[arr,:] .= div[arr,:] .+
    reshape(vC1[:,:,i,:]
    ,OP*OP*size(i,1),nz);
end
div=div./CG.M;
divCG=reshape(div[reshape(CG.Glob,OP*OP*NF,1),:]
  ,OP,OP,NF,nz);
return divCG
end
