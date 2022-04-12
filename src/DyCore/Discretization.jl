mutable struct CGStruct
    OrdPoly::Int
    OrdPolyZ::Int
    Glob::Array{Int, 2}
    Stencil::Array{Float64, 2}
    NumG::Int
    NumI::Int
    w::Array{Float64, 1}
    xw::Array{Float64, 1}
    xwZ::Array{Float64, 1}
    DW::Array{Float64, 2}
    DWT::Array{Float64, 2}
    DS::Array{Float64, 2}
    DST::Array{Float64, 2}
    M::Array{Float64, 2}
    MW::Array{Float64, 2}
end
function CGStruct()
 OrdPoly=0
OrdPolyZ=0
Glob=zeros(0,0)
Stencil=zeros(0,0)
NumG=0
NumI=0
w=zeros(0)
xw=zeros(0)
xwZ=zeros(0)
DW=zeros(0,0)
DWT=zeros(0,0)
DS=zeros(0,0)
DST=zeros(0,0)
M=zeros(0,0)
MW=zeros(0,0)
 return CGStruct(
    OrdPoly,
    OrdPolyZ,
    Glob,
    Stencil,
    NumG,
    NumI,
    w,
    xw,
    xwZ,
    DW,
    DWT,
    DS,
    DST,
    M,
    MW
 )
end 

function Discretization(OrdPoly,OrdPolyZ,Jacobi,Global)
# Discretization
OP=OrdPoly+1
OPZ=OrdPolyZ+1
nz=Global.Grid.nz;
NF=Global.Grid.NumFaces

CG = CGStruct()
CG.OrdPoly=OrdPoly;
CG.OrdPolyZ=OrdPolyZ;

(CG.w,CG.xw)=GaussLobattoQuad(CG.OrdPoly);
(wZ,CG.xwZ)=GaussLobattoQuad(CG.OrdPolyZ);
(CG.DW,CG.DS)=DerivativeMatrixSingle(CG.OrdPoly);
CG.DST=CG.DS'
CG.DWT=CG.DW'
(CG.Glob,CG.NumG,CG.NumI,CG.Stencil) =
  NumberingFemCG(Global.Grid,OrdPoly);


dXdxIF = Global.Metric.dXdxIF
dXdxIC = Global.Metric.dXdxIC
J = Global.Metric.J;
JC = Global.Metric.JC;
JF = Global.Metric.JF;
lat = Global.Metric.lat;
X = Global.Metric.X;
dXdx   = zeros(OP,OP,OPZ,nz,3,3,NF)
dXdxI  = zeros(OP,OP,OPZ,nz,3,3,NF)

for iF=1:Global.Grid.NumFaces
  for iz=1:nz
    zI=[Global.Grid.z[iz],Global.Grid.z[iz+1]];
    (X_Fz,J_Fz,dXdx_Fz,dXdxI_Fz,lat_Fz)=Jacobi(CG,Global.Grid.Faces[iF],zI,Topo,Global.Grid.Topography);
    X[:,:,:,:,iz,iF]=X_Fz;
    J[:,:,:,iz,iF]=J_Fz;
    dXdx[:,:,:,iz,:,:,iF]=reshape(dXdx_Fz,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
    dXdxI[:,:,:,iz,:,:,iF]=reshape(dXdxI_Fz,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
    lat[:,:,iF]=lat_Fz;
  end
end
@views @. JC=0.5*(J[:,:,2,:,:] + J[:,:,1,:,:])
@views @. JF[:,:,1,:] = J[:,:,1,1,:]
@views @. JF[:,:,2:nz,:] = 0.5*(J[:,:,2,1:nz-1,:] + J[:,:,1,2:nz,:])
@views @. JF[:,:,nz+1,:] = J[:,:,2,nz,:]

AverageFB!(JF,J);
dXdxIC .= 0;
dXdxIF .= 0;
for i=1:3
  for j=1:3
    @views @. dXdxIC[:,:,:,i,j,:] = 0.5*(dXdxI[:,:,1,:,i,j,:] + dXdxI[:,:,2,:,i,j,:])
    @views @. dXdxIF[:,:,1,i,j,:] = dXdxI[:,:,1,1,i,j,:]
    @views @. dXdxIF[:,:,2:nz,i,j,:] = 0.5*(dXdxI[:,:,2,1:nz-1,i,j,:] + dXdxI[:,:,1,2:nz,i,j,:])
    @views @. dXdxIF[:,:,nz+1,i,j,:] = dXdxI[:,:,2,nz,i,j,:]
  end
end

(CG.M,CG.MW)=MassCG(CG,Global);
Global.latN=zeros(CG.NumG);
lat = Global.Metric.lat
latN = Global.latN
OP=CG.OrdPoly+1;
for iF=1:NF
  iG=0
  for iP=1:OP
    for jP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
      latN[ind] = latN[ind] + lat[iP,jP,iF]
    end
  end
end


Global.latN=Global.latN./CG.M[1,:];
return (CG,Global)
end
