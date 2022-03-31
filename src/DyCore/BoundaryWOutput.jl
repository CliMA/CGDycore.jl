function BoundaryWOutput!(W,U,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
dXdxIC = Param.cache.dXdxIC;
JC = Param.cache.JC;

v1CG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),1,Param.uPos],OP,OP,NF);
v2CG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),1,Param.vPos],OP,OP,NF,1);

wCG = Param.wCG[:,:,:,1]
@views wCG .= -(dXdxIC[:,:,:,1,3,1] .* v1CG[:,:,:,1] .+
   dXdxIC[:,:,:,1,3,2] .* v2CG[:,:,:,1]) ./ dXdxIC[:,:,:,1,3,3] .* JC[:,:,:,1];
W .= 0
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  @views arr = reshape(CG.Glob[:,i,:],OP*OP*size(i,1))
  @views mat = reshape(wCG[:,:,i] ,OP*OP*size(i,1),1)
  @views W[arr,:,:] = W[arr,:,:] .+ mat;
end
W=W./CG.M[:,1];
end
