function FGrad3Vec!(F,cCG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;

D1cCG = @view Param.CacheC[:,:,:,:,1]
D2cCG = @view Param.CacheC[:,:,:,:,2]
D3cCG = @view Param.CacheC[:,:,:,:,3]
D3cCGE = @view Param.CacheC[:,:,:,:,4]
mul!(reshape(D1cCG,OP,OP*NF*nz),CG.DS,reshape(cCG,OP,OP*nz*NF))
mul!(reshape(PermutedDimsArray(D2cCG,(2,1,3,4)),OP,OP*NF*nz),CG.DS,reshape(PermutedDimsArray(cCG,(2,1,3,4)),OP,OP*nz*NF))
@views D3cCG[:,:,:,1:end-1] .= 0.5*(cCG[:,:,:,2:end] .- cCG[:,:,:,1:end-1])

@views F[:,:,:,1:end-1,Param.wPos] .= F[:,:,:,1:end-1,Param.wPos] .- Param.dXdxIF[:,:,:,2:end-1,3,3].*D3cCG[:,:,:,1:end-1];
if nz>1
@views  D3cCGE[:,:,:,1]=D3cCG[:,:,:,1];
@views  D3cCGE[:,:,:,2:end-1] .= 0.5*(D3cCG[:,:,:,1:end-2] .+ D3cCG[:,:,:,2:end-1]);
@views  D3cCGE[:,:,:,end]=D3cCG[:,:,:,end-1];
else
    @views  D3cCGE[:,:,:,1] .= 0
end

@views F[:,:,:,:,Param.uPos] .= F[:,:,:,:,Param.uPos] .- (Param.dXdxIC[:,:,:,:,1,1].*D1cCG .+
  Param.dXdxIC[:,:,:,:,2,1].*D2cCG .+
  Param.dXdxIC[:,:,:,:,3,1].*D3cCGE)  

@views F[:,:,:,:,Param.vPos] .= F[:,:,:,:,Param.vPos] .- (Param.dXdxIC[:,:,:,:,1,2].*D1cCG .+
  Param.dXdxIC[:,:,:,:,2,2].*D2cCG .+
  Param.dXdxIC[:,:,:,:,3,2].*D3cCGE) 

end


function FGrad3RhoVec!(F,cCG,RhoCG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;

D1cCG = @view Param.CacheC[:,:,:,:,1]
D2cCG = @view Param.CacheC[:,:,:,:,2]
D3cCG = @view Param.CacheC[:,:,:,:,3]
D3cCGE = @view Param.CacheC[:,:,:,:,4]

mul!(reshape(D1cCG,OP,OP*NF*nz),CG.DS,reshape(cCG,OP,OP*nz*NF))
mul!(reshape(PermutedDimsArray(D2cCG,(2,1,3,4)),OP,OP*NF*nz),CG.DS,reshape(PermutedDimsArray(cCG,(2,1,3,4)),OP,OP*nz*NF))
@views D3cCG[:,:,:,1:end-1] .= 0.5*(cCG[:,:,:,2:end] .- cCG[:,:,:,1:end-1])

@views F[:,:,:,1:end-1,Param.wPos] .= F[:,:,:,1:end-1,Param.wPos] .- 
    Param.dXdxIF[:,:,:,2:end-1,3,3].*D3cCG[:,:,:,1:end-1] ./ (0.5*(RhoCG[:,:,:,1:end-1]+RhoCG[:,:,:,2:end]))
if nz>1
@views  D3cCGE[:,:,:,1] .= D3cCG[:,:,:,1];
@views  D3cCGE[:,:,:,2:end-1] .= 0.5*(D3cCG[:,:,:,1:end-2] .+ D3cCG[:,:,:,2:end-1]);
@views  D3cCGE[:,:,:,end] .= D3cCG[:,:,:,end-1];
else
    @views  D3cCGE[:,:,:,1] .= 0
end

@views F[:,:,:,:,Param.uPos] .= F[:,:,:,:,Param.uPos] .- (Param.dXdxIC[:,:,:,:,1,1].*D1cCG .+
  Param.dXdxIC[:,:,:,:,2,1].*D2cCG .+
  Param.dXdxIC[:,:,:,:,3,1].*D3cCGE) ./ RhoCG  

@views F[:,:,:,:,Param.vPos] .= F[:,:,:,:,Param.vPos] .- (Param.dXdxIC[:,:,:,:,1,2].*D1cCG .+
  Param.dXdxIC[:,:,:,:,2,2].*D2cCG .+
  Param.dXdxIC[:,:,:,:,3,2].*D3cCGE) ./ RhoCG

end


