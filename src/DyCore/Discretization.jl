Base.@kwdef mutable struct CGStruct
    OrdPoly = nothing
    OrdPolyZ = nothing
    Faces = nothing
    NumG = nothing
    NumI = nothing
    Glob = nothing
    FaceGlob = nothing
    Stencil = nothing
    w = nothing
    xw = nothing
    wX = nothing
    xwX = nothing
    wY = nothing
    xwY = nothing
    wZ = nothing
    xwZ = nothing
    DW = nothing
    DS = nothing
    DWZ = nothing
    DSZ = nothing
    M = nothing
    MW = nothing
end

function Discretization(OrdPoly,OrdPolyZ,Jacobi,Param)
# Discretization

CG = CGStruct(;)
CG.OrdPoly=OrdPoly;
CG.OrdPolyZ=OrdPolyZ;
nz=Param.Grid.nz;
dXdxIF = Param.cache.dXdxIF
dXdxIC = Param.cache.dXdxIC

(CG.Faces,CG.NumG,CG.NumI,CG.Glob,CG.FaceGlob,CG.Stencil) =
  NumberingFemCG(Param.Grid,OrdPoly);
Param.J=zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,Param.Grid.NumFaces,nz);
Param.lat=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces);
Param.X=zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,3,Param.Grid.NumFaces,nz);
Param.dXdx=zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,Param.Grid.NumFaces,nz,3,3);
Param.dXdxI=zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,Param.Grid.NumFaces,nz,3,3);
(CG.w,CG.xw)=GaussLobattoQuad(CG.OrdPoly);
(CG.wX,CG.xwX)=GaussLobattoQuad(CG.OrdPoly);
(CG.wY,CG.xwY)=GaussLobattoQuad(CG.OrdPoly);
(CG.wZ,CG.xwZ)=GaussLobattoQuad(CG.OrdPolyZ);
(CG.DW,CG.DS)=DerivativeMatrixSingle(CG.OrdPoly);
(CG.DWZ,CG.DSZ)=DerivativeMatrixSingle(CG.OrdPolyZ);


for iF=1:Param.Grid.NumFaces
  for iz=1:nz
    zI=[Param.Grid.z[iz],Param.Grid.z[iz+1]];
    (X,J,dXdx,dXdxI,lat)=Jacobi(CG,Param.Grid.Faces[iF],zI,Topo,Param);
    Param.X[:,:,:,:,iF,iz]=X;
    Param.J[:,:,:,iF,iz]=J;
    Param.dXdx[:,:,:,iF,iz,:,:]=reshape(dXdx,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
    Param.dXdxI[:,:,:,iF,iz,:,:]=reshape(dXdxI,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
    Param.lat[:,:,iF]=lat;
  end
end
Param.JC=Average(Param.J);
Param.JF=AverageFB(Param.J);
dXdxIC .= 0;
dXdxIF .= 0;
for i=1:3
  for j=1:3
    dXdxIC[:,:,:,:,i,j] .= Average(Param.dXdxI[:,:,:,:,:,i,j]);
    dXdxIF[:,:,:,:,i,j] .= AverageFB(Param.dXdxI[:,:,:,:,:,i,j]);
  end
end

(CG.M,CG.MW)=MassCG(CG,Param);
Param.latN=zeros(CG.NumG,1);
OP=CG.OrdPoly+1;
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:,:] ,OP*OP*size(i,1),1)
  mat = Param.lat[:,:,i] .*Param.JC[:,:,i,1]
  Param.latN[arr,:,:]= Param.latN[arr,:,:] + reshape(mat ,OP*OP*size(i,1),1);
end
Param.latN=Param.latN./CG.M[:,1];
return (CG,Param)
end
