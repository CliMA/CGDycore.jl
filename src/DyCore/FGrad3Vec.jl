function FGrad3Vec!(F,cCG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;

D1cCG = Param.CacheC1
D2cCG = Param.CacheC2
D3cCG = Param.CacheC3
D3cCGE =Param.CacheC4
mul!(reshape(D1cCG,OP,OP*NF*nz),CG.DS,reshape(cCG,OP,OP*nz*NF))
mul!(reshape(PermutedDimsArray(D2cCG,(2,1,3,4)),OP,OP*NF*nz),CG.DS,reshape(PermutedDimsArray(cCG,(2,1,3,4)),OP,OP*nz*NF))
@views D3cCG[:,:,:,1:nz-1] .= 0.5*(cCG[:,:,:,2:nz] .- cCG[:,:,:,1:nz-1])

@views F[:,:,:,1:nz-1,Param.wPos] .= F[:,:,:,1:nz-1,Param.wPos] .- Param.dXdxIF33.*D3cCG[:,:,:,1:nz-1];

if nz>1
@views  D3cCGE[:,:,:,1] .= D3cCG[:,:,:,1];
@views  D3cCGE[:,:,:,2:end-1] .= 0.5.*(D3cCG[:,:,:,1:end-2] .+ D3cCG[:,:,:,2:end-1]);
@views  D3cCGE[:,:,:,end] .= D3cCG[:,:,:,end-1];
else
    @views  D3cCGE[:,:,:,1] .= 0
end

@views F[:,:,:,:,Param.uPos] .= F[:,:,:,:,Param.uPos] .- (Param.dXdxIC11.*D1cCG .+
  Param.dXdxIC21.*D2cCG .+
  Param.dXdxIC31.*D3cCGE)

@views F[:,:,:,:,Param.vPos] .= F[:,:,:,:,Param.vPos] .- (Param.dXdxIC12.*D1cCG .+
  Param.dXdxIC22.*D2cCG .+ Param.dXdxIC32.*D3cCGE)

end


function FGrad3RhoVec!(F,cCG,RhoCG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;

D1cCG = Param.CacheC1
D2cCG = Param.CacheC2
D3cCG = Param.CacheC3
D3cCGE = Param.CacheC4

mul!(reshape(D1cCG,OP,OP*NF*nz),CG.DS,reshape(cCG,OP,OP*nz*NF))
mul!(reshape(PermutedDimsArray(D2cCG,(2,1,3,4)),OP,OP*NF*nz),CG.DS,reshape(PermutedDimsArray(cCG,(2,1,3,4)),OP,OP*nz*NF))
@views D3cCG[:,:,:,1:end-1] .= 0.5.*(cCG[:,:,:,2:end] .- cCG[:,:,:,1:end-1])

@views F[:,:,:,1:nz-1,Param.wPos] .= F[:,:,:,1:nz-1,Param.wPos] .-
    Param.dXdxIF33.*D3cCG[:,:,:,1:nz-1] ./ (0.5.*(RhoCG[:,:,:,1:nz-1]+RhoCG[:,:,:,2:nz]))
if nz>1
@views  D3cCGE[:,:,:,1] .= D3cCG[:,:,:,1];
@views  D3cCGE[:,:,:,2:end-1] .= 0.5.*(D3cCG[:,:,:,1:nz-2] .+ D3cCG[:,:,:,2:nz-1]);
@views  D3cCGE[:,:,:,end] .= D3cCG[:,:,:,nz-1];
else
    @views  D3cCGE[:,:,:,1] .= 0
end

@views F[:,:,:,:,Param.uPos] .= F[:,:,:,:,Param.uPos] .- (Param.dXdxIC11.*D1cCG .+
  Param.dXdxIC21.*D2cCG .+
  Param.dXdxIC31.*D3cCGE) ./ RhoCG

@views F[:,:,:,:,Param.vPos] .= F[:,:,:,:,Param.vPos] .- (Param.dXdxIC12.*D1cCG .+
  Param.dXdxIC22.*D2cCG .+
  Param.dXdxIC32.*D3cCGE) ./ RhoCG

end

function FGrad3Vec(cCG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
D1cCG=reshape(
  CG.DS*reshape(cCG,OP,OP*NF*nz)
  ,OP,OP,NF,nz);
D2cCG= permute(
  reshape(
  CG.DS *reshape(permute(cCG
  ,[2 1 3 4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2 1 3 4]);

D3cCG=0.5*(cCG[:,:,:,2:nz]-cCG[:,:,:,1:nz-1]);

gradCG=zeros(OP,OP,NF,nz,3);
gradCG[:,:,:,1:nz-1,3]=Param.dXdxIF[:,:,:,2:nz,3,3].*D3cCG;
D3cCGE=zeros(OP,OP,NF,nz);
if nz>1
  D3cCGE[:,:,:,1]=D3cCG[:,:,:,1];
  D3cCGE[:,:,:,2:nz-1]=0.5*(D3cCG[:,:,:,1:end-1]+D3cCG[:,:,:,2:end]);
  D3cCGE[:,:,:,nz]=D3cCG[:,:,:,nz-1];
end

gradCG[:,:,:,:,1]=Param.dXdxIC[:,:,:,:,1,1].*D1cCG+
  Param.dXdxIC[:,:,:,:,2,1].*D2cCG+
  Param.dXdxIC[:,:,:,:,3,1].*D3cCGE;
gradCG[:,:,:,:,2]=Param.dXdxIC[:,:,:,:,1,2].*D1cCG+
  Param.dXdxIC[:,:,:,:,2,2].*D2cCG+
  Param.dXdxIC[:,:,:,:,3,2].*D3cCGE;

return gradCG
end

