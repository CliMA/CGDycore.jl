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
J = Param.cache.J;
JC = Param.cache.JC;
lat = Param.cache.lat;
X = Param.cache.X;
dXdx = Param.cache.dXdx;
dXdxI = Param.cache.dXdxI;

J .= 0;
lat .= 0;
X .= 0;
dXdx .= 0;
dXdxI .= 0;
(CG.w,CG.xw)=GaussLobattoQuad(CG.OrdPoly);
(CG.wX,CG.xwX)=GaussLobattoQuad(CG.OrdPoly);
(CG.wY,CG.xwY)=GaussLobattoQuad(CG.OrdPoly);
(CG.wZ,CG.xwZ)=GaussLobattoQuad(CG.OrdPolyZ);
(CG.DW,CG.DS)=DerivativeMatrixSingle(CG.OrdPoly);
(CG.DWZ,CG.DSZ)=DerivativeMatrixSingle(CG.OrdPolyZ);


for iF=1:Param.Grid.NumFaces
  for iz=1:nz
    zI=[Param.Grid.z[iz],Param.Grid.z[iz+1]];
    (X_Fz,J_Fz,dXdx_Fz,dXdxI_Fz,lat_Fz)=Jacobi(CG,Param.Grid.Faces[iF],zI,Topo,Param);
    X[:,:,:,:,iF,iz]=X_Fz;
    J[:,:,:,iF,iz]=J_Fz;
    dXdx[:,:,:,iF,iz,:,:]=reshape(dXdx_Fz,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
    dXdxI[:,:,:,iF,iz,:,:]=reshape(dXdxI_Fz,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
    lat[:,:,iF]=lat_Fz;
  end
end
JC .= Average(J);
Param.JF=AverageFB(J);
dXdxIC .= 0;
dXdxIF .= 0;
for i=1:3
  for j=1:3
    dXdxIC[:,:,:,:,i,j] .= Average(dXdxI[:,:,:,:,:,i,j]);
    dXdxIF[:,:,:,:,i,j] .= AverageFB(dXdxI[:,:,:,:,:,i,j]);
  end
end
Param.dXdxIC11 = view(Param.dXdxIC,:,:,:,:,1,1)
Param.dXdxIC12 = view(Param.dXdxIC,:,:,:,:,1,2)
Param.dXdxIC21 = view(Param.dXdxIC,:,:,:,:,2,1)
Param.dXdxIC22 = view(Param.dXdxIC,:,:,:,:,2,2)
Param.dXdxIC31 = view(Param.dXdxIC,:,:,:,:,3,1)
Param.dXdxIC32 = view(Param.dXdxIC,:,:,:,:,3,2)
Param.dXdxIC33 = view(Param.dXdxIC,:,:,:,:,3,3)
Param.dXdxIF11 = view(Param.dXdxIF,:,:,:,:,1,1)
Param.dXdxIF12 = view(Param.dXdxIF,:,:,:,:,1,2)
Param.dXdxIF21 = view(Param.dXdxIF,:,:,:,:,2,1)
Param.dXdxIF22 = view(Param.dXdxIF,:,:,:,:,2,2)
Param.dXdxIF31 = view(Param.dXdxIF,:,:,:,:,3,1)
Param.dXdxIF32 = view(Param.dXdxIF,:,:,:,:,3,2)
Param.dXdxIC33 = view(Param.dXdxIC,:,:,:,:,3,3)

(CG.M,CG.MW)=MassCG(CG,Param);
Param.latN=zeros(CG.NumG,1);
lat = Param.cache.lat
OP=CG.OrdPoly+1;
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:,:] ,OP*OP*size(i,1),1)
  mat = lat[:,:,i] .*JC[:,:,i,1]
  Param.latN[arr,:,:]= Param.latN[arr,:,:] + reshape(mat ,OP*OP*size(i,1),1);
end
Param.latN=Param.latN./CG.M[:,1];
return (CG,Param)
end
