function FDiv3Vec(cCG,v1CG,v2CG,v3CG,CG,Param)
# noch uberarbeiten
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
dXdxIF = Param.cache.dXdxIF;
dXdxIC = Param.cache.dXdxIC;
# Contravariant components
v1Con=v1CG.* dXdxIC[:,:,:,:,1,1] +
   v2CG.* dXdxIC[:,:,:,:,1,2];
v2Con=v1CG.* dXdxIC[:,:,:,:,2,1] +
   v2CG.* dXdxIC[:,:,:,:,2,2];


v3Con=0.5*((v1CG[:,:,:,1:end-1]+v1CG[:,:,:,2:end]).*
        dXdxIF[:,:,:,2:nz,3,1]+
      (v2CG[:,:,:,1:end-1]+v2CG[:,:,:,2:end]).*
        dXdxIF[:,:,:,2:nz,3,2]) +
       v3CG[:,:,:,2:nz].*
        dXdxIF[:,:,:,2:nz,3,3];

Dv1Con=reshape(
  CG.DS*reshape(v1Con.*cCG,OP,OP*NF*nz),
   OP,OP,NF,nz);
Dv2Con=permute(
  reshape(
  CG.DS*reshape(
  permute(v2Con.*cCG,
  [2,1,3,4]),
  OP,OP*NF*nz),
  OP,OP,NF,nz),
  [2,1,3,4]);
v3Con=0.5*(cCG[:,:,:,1:nz-1]+cCG[:,:,:,2:nz]).*v3Con;
Dv3Con=zeros(OP,OP,NF,nz);
if nz>1
  Dv3Con[:,:,:,1]=0.5*v3Con[:,:,:,1];
  Dv3Con[:,:,:,2:nz-1]=0.5*(v3Con[:,:,:,2:nz-1]-v3Con[:,:,:,1:nz-2]);
  Dv3Con[:,:,:,nz]=-0.5*v3Con[:,:,:,nz-1]; # 0.5 Metric
end
divCG=Dv1Con+Dv2Con+Dv3Con;

return divCG
end
