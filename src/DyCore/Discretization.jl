mutable struct CGStruct
    OrdPoly::Int
    OrdPolyZ::Int
    Glob::Array{Int, 3}
    Stencil::Array{Int, 2}
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
    MMass::Array{Float64, 2}
    MW::Array{Float64, 2}
    Boundary::Array{Int, 1}
end
function CGStruct()
 OrdPoly=0
OrdPolyZ=0
Glob=zeros(0,0,0)
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
MMass=zeros(0,0)
MW=zeros(0,0)
Boundary=zeros(0)
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
    MMass,
    MW,
    Boundary,
 )
end 

function Discretization(OrdPoly,OrdPolyZ,Jacobi,Global)
# Discretization
  Grid = Global.Grid
  OP=OrdPoly+1
  OPZ=OrdPolyZ+1
  nz=Grid.nz;
  NF=Grid.NumFaces

  CG = CGStruct()
  CG.OrdPoly=OrdPoly;
  CG.OrdPolyZ=OrdPolyZ;

  (CG.w,CG.xw)=GaussLobattoQuad(CG.OrdPoly);
  (wZ,CG.xwZ)=GaussLobattoQuad(CG.OrdPolyZ);
  (CG.DW,CG.DS)=DerivativeMatrixSingle(CG.OrdPoly);
  CG.DST=CG.DS'
  CG.DWT=CG.DW'
  (CG.Glob,CG.NumG,CG.NumI,CG.Stencil) =
    NumberingFemCG(Grid,OrdPoly);


  dXdxIF = Global.Metric.dXdxIF
  dXdxIC = Global.Metric.dXdxIC
  nS = Global.Metric.nS
  dz = Global.Metric.dz
  J = Global.Metric.J;
  JC = Global.Metric.JC;
  JF = Global.Metric.JF;
  lat = Global.Metric.lat;
  X = Global.Metric.X;
  dXdx   = zeros(OP,OP,OPZ,nz,3,3,NF)
  dXdxI  = zeros(OP,OP,OPZ,nz,3,3,NF)

  for iF=1:Grid.NumFaces
    for iz=1:nz
      zI=[Grid.z[iz],Grid.z[iz+1]];
      (X_Fz,J_Fz,dXdx_Fz,dXdxI_Fz,lat_Fz)=Jacobi(CG,Grid.Faces[iF],zI,Topo,Grid.Topography);
      X[:,:,:,:,iz,iF]=X_Fz;
      J[:,:,:,iz,iF]=J_Fz;
      dXdx[:,:,:,iz,:,:,iF]=reshape(dXdx_Fz,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
      dXdxI[:,:,:,iz,:,:,iF]=reshape(dXdxI_Fz,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
      lat[:,:,iF]=lat_Fz;
    end
#   Surface normal
    @views @. nS[:,:,1,iF] = dXdxI[:,:,1,1,3,1,iF] / sqrt(dXdxI[:,:,1,1,3,1,iF] * dXdxI[:,:,1,1,3,1,iF] +
      dXdxI[:,:,1,1,3,2,iF] * dXdxI[:,:,1,1,3,2,iF] + dXdxI[:,:,1,1,3,3,iF] * dXdxI[:,:,1,1,3,3,iF])
    @views @. nS[:,:,2,iF] = dXdxI[:,:,1,1,3,2,iF] / sqrt(dXdxI[:,:,1,1,3,1,iF] * dXdxI[:,:,1,1,3,1,iF] +
      dXdxI[:,:,1,1,3,2,iF] * dXdxI[:,:,1,1,3,2,iF] + dXdxI[:,:,1,1,3,3,iF] * dXdxI[:,:,1,1,3,3,iF])
    @views @. nS[:,:,3,iF] = dXdxI[:,:,1,1,3,3,iF] / sqrt(dXdxI[:,:,1,1,3,1,iF] * dXdxI[:,:,1,1,3,1,iF] +
      dXdxI[:,:,1,1,3,2,iF] * dXdxI[:,:,1,1,3,2,iF] + dXdxI[:,:,1,1,3,3,iF] * dXdxI[:,:,1,1,3,3,iF])
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

  @views @. dz = 2.0 * JC / dXdxIC[:,:,:,3,3,:]



  (CG.M,CG.MW,CG.MMass)=MassCG(CG,Global);
  Global.latN=zeros(CG.NumG);
  lat = Global.Metric.lat
  latN = Global.latN
  OP=CG.OrdPoly+1;
  for iF=1:NF
    for jP=1:OP
      for iP=1:OP
        ind=CG.Glob[iP,jP,iF]
        latN[ind] = latN[ind] + lat[iP,jP,iF]
      end
    end
  end
  Global.latN=Global.latN./CG.M[1,:];
# Boundary nodes
  for iF = 1 : NF
    Side = 0
    for iE in Grid.Faces[iF].E
       Side += 1 
       if Grid.Edges[iE].F[1] == 0 || Grid.Edges[iE].F[2] == 0
         if Side == 1
           for i in CG.Glob[1:OP-1,1,iF]   
             push!(CG.Boundary,i)
           end
         elseif Side == 2  
           for i in CG.Glob[OP,1:OP-1,iF]
             push!(CG.Boundary,i)
           end
         elseif Side == 3  
           for i in CG.Glob[2:OP,OP,iF]
             push!(CG.Boundary,i)
           end
         elseif Side == 4  
           for i in CG.Glob[1,2:OP,iF]
             push!(CG.Boundary,i)
           end
         end  
       end  
    end
  end  
  return (CG,Global)
end
