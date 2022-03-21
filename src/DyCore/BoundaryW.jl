BoundaryW(v1CG,v2CG,CG,Param) =
  BoundaryW(v1CG,v2CG,CG,Param.cache.dXdxIC)

function BoundaryW(v1CG,v2CG,CG,dXdxIC::AbstractArray)
wCG = -(dXdxIC[:,:,:,1,3,1].*v1CG[:,:,:,1]+
  dXdxIC[:,:,:,1,3,2].*v2CG[:,:,:,1])./
  dXdxIC[:,:,:,1,3,3];
  return wCG
end