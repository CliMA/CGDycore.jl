function BoundaryWOutput!(W,U,CG,Global)
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
dXdxIC = Global.Metric.dXdxIC;
JC = Global.Metric.JC;
@views v1CG = Global.Cache.v1CG[:,:,1]
@views v2CG = Global.Cache.v2CG[:,:,1]
@views wCG = Global.Cache.wCG[:,:,1]

W .= 0
@inbounds for iF=1:NF
  iG=0
  @inbounds for iP=1:OP
    @inbounds for jP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
      v1CG[iP,jP] = U[iG,1,Global.Model.uPos]
      v2CG[iP,jP] = U[iG,1,Global.Model.vPos]
    end
  end

  @views @. wCG = -(dXdxIC[:,:,1,3,1,iF] * v1CG[:,:] +
    dXdxIC[:,:,1,3,2,iF] .* v2CG[:,:,]) / dXdxIC[:,:,1,3,3,iF] * JC[:,:,1,iF];
    
  iG=0
  for iP=1:OP
    for jP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
      W[ind] = W[ind] + wCG[iP,jP]
    end
  end
end
W .= W./CG.M[1,:];
end
