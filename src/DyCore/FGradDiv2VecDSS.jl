function FGradDiv2VecDSS!(grad1CG,grad2CG,v1CG,v2CG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;

vCon = Param.CacheC1
DvCon = Param.CacheC2
vC1 = Param.CacheC3
D1vC1 = vCon
D2vC1 = DvCon
grad1TempCG = Param.CacheC3
grad2TempCG = Param.CacheC4
grad1 = Param.Cache1
grad2 = Param.Cache2

vCon .= v1CG.*Param.dXdxIC11 .+ v2CG.*Param.dXdxIC12 ;
mul!(reshape(vC1,OP,OP*NF*nz),CG.DS,reshape(vCon,OP,OP*nz*NF))
vCon .= v1CG.*Param.dXdxIC21 .+ v2CG.*Param.dXdxIC22;
mul!(reshape(PermutedDimsArray(DvCon,(2,1,3,4)),OP,OP*NF*nz),CG.DS,reshape(PermutedDimsArray(vCon,(2,1,3,4)),OP,OP*nz*NF))
vC1 .= vC1 .+ DvCon

mul!(reshape(D1vC1,OP,OP*NF*nz),CG.DW,reshape(vC1,OP,OP*nz*NF))
mul!(reshape(PermutedDimsArray(D2vC1,(2,1,3,4)),OP,OP*NF*nz),CG.DW,reshape(PermutedDimsArray(vC1,(2,1,3,4)),OP,OP*nz*NF))


grad1TempCG .= (Param.dXdxIC11.*D1vC1 .+ Param.dXdxIC21.*D2vC1)./Param.JC;
grad2TempCG .= (Param.dXdxIC12.*D1vC1 .+ Param.dXdxIC22.*D2vC1)./Param.JC;

grad1 .= 0
grad2 .= 0
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:],OP*OP*size(i,1))
  grad1[arr,:,:] .= grad1[arr,:,:] .+
    reshape(grad1TempCG[:,:,i,:]
    ,OP*OP*size(i,1),nz);
  grad2[arr,:,:] .= grad2[arr,:,:] .+
    reshape(grad2TempCG[:,:,i,:]
    ,OP*OP*size(i,1),nz);
end

grad1.=grad1./CG.M;
grad2.=grad2./CG.M;
grad1CG .= reshape(grad1[reshape(CG.Glob,OP*OP*NF,1),:]
  ,OP,OP,NF,nz);
grad2CG .= reshape(grad2[reshape(CG.Glob,OP*OP*NF,1),:]
  ,OP,OP,NF,nz);

end

function FGradDiv2VecDSS(v1CG,v2CG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
vC1=reshape(
  CG.DS*reshape(v1CG.*Param.dXdxIC[:,:,:,:,1,1] +
  v2CG.*Param.dXdxIC[:,:,:,:,1,2]
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz) +
  permute(
  reshape(
  CG.DS*reshape(
  permute(v1CG.*Param.dXdxIC[:,:,:,:,2,1] +
  v2CG.*Param.dXdxIC[:,:,:,:,2,2]
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
gradCG[:,:,:,:,1]=(Param.dXdxIC[:,:,:,:,1,1].*D1cCG +
  Param.dXdxIC[:,:,:,:,2,1].*D2cCG)./Param.JC;
gradCG[:,:,:,:,2]=(Param.dXdxIC[:,:,:,:,1,2].*D1cCG +
  Param.dXdxIC[:,:,:,:,2,2].*D2cCG)./Param.JC;


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

