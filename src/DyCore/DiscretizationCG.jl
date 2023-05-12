mutable struct CGStruct
    OrdPoly::Int
    OrdPolyZ::Int
    Glob::Array{Int, 3}
    Stencil::Array{Int, 2}
    NumG::Int
    NumI::Int
    w::Array{Float64, 1}
    xw::Array{Float64, 1}
    xe::Array{Float64, 1}
    IntXE2F::Array{Float64, 2}
    xwZ::Array{Float64, 1}
    IntZE2F::Array{Float64, 2}
    DW::Array{Float64, 2}
    DWT::Array{Float64, 2}
    DS::Array{Float64, 2}
    DST::Array{Float64, 2}
    DSZ::Array{Float64, 2}
    S::Array{Float64, 2}
    M::Array{Float64, 2}
    MMass::Array{Float64, 2}
    MW::Array{Float64, 2}
    Boundary::Array{Int, 1}
    MasterSlave::Array{Int, 1}
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
xe=zeros(0)
IntXE2F=zeros(0,0)
xwZ=zeros(0)
IntZE2F=zeros(0,0)
DW=zeros(0,0)
DWT=zeros(0,0)
DS=zeros(0,0)
DST=zeros(0,0)
DSZ=zeros(0,0)
S=zeros(0,0)
M=zeros(0,0)
MMass=zeros(0,0)
MW=zeros(0,0)
Boundary=zeros(0)
MasterSlave = zeros(0)
 return CGStruct(
    OrdPoly,
    OrdPolyZ,
    Glob,
    Stencil,
    NumG,
    NumI,
    w,
    xw,
    xe,
    IntXE2F,
    xwZ,
    IntZE2F,
    DW,
    DWT,
    DS,
    DST,
    DSZ,
    S,
    M,
    MMass,
    MW,
    Boundary,
    MasterSlave,
 )
end 

function DiscretizationCG(OrdPoly,OrdPolyZ,Jacobi,Global)
  DiscretizationCG(OrdPoly,OrdPolyZ,Jacobi,Global,zeros(OrdPoly+1,OrdPoly+1,Global.Grid.NumFaces))
end  

function DiscretizationCG(OrdPoly,OrdPolyZ,Jacobi,Global,zs)
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
  CG.xe = zeros(OrdPoly+1)
  CG.xe[1] = -1.0
  for i = 2 : OrdPoly
    CG.xe[i] = CG.xe[i-1] + 2.0/OrdPoly
  end
  CG.xe[OrdPoly+1] = 1.0

  CG.IntXE2F = zeros(OrdPoly+1,OrdPoly+1)
  for j = 1 : OrdPoly + 1
    for i = 1 : OrdPoly +1
      CG.IntXE2F[i,j] = Lagrange(CG.xw[i],CG.xe,j)
    end
  end

  CG.IntZE2F = zeros(OrdPolyZ+1,OrdPolyZ+1)
  for j = 1 : OrdPolyZ + 1
    for i = 1 : OrdPolyZ +1
      CG.IntZE2F[i,j] = Lagrange(CG.xwZ[i],CG.xwZ,j)
    end
  end

  (CG.DW,CG.DS)=DerivativeMatrixSingle(CG.OrdPoly);
  CG.DST=CG.DS'
  CG.DWT=CG.DW'

  Q = diagm(CG.w) * CG.DS
  CG.S = Q - Q'
  (DWZ,CG.DSZ)=DerivativeMatrixSingle(CG.OrdPolyZ);
  (CG.Glob,CG.NumG,CG.NumI,CG.Stencil,CG.MasterSlave) =
    NumberingFemCG(Grid,OrdPoly);  


  dXdxI = Global.Metric.dXdxI
  nS = Global.Metric.nS
  FS = Global.Metric.FS
  Global.Metric.dz = zeros(nz,CG.NumG)
  Global.Metric.zP = zeros(nz,CG.NumG)
  dz = Global.Metric.dz
  zP = Global.Metric.zP
  J = Global.Metric.J;
  lat = Global.Metric.lat;
  X = Global.Metric.X;
  dXdx   = zeros(OP,OP,OPZ,nz,3,3,NF)
  dXdxILoc  = zeros(OP,OP,OPZ,nz,3,3,NF)

  for iF=1:Grid.NumFaces
    for iz=1:nz
      zI=[Grid.z[iz],Grid.z[iz+1]];
      (X_Fz,J_Fz,dXdx_Fz,dXdxI_Fz,lat_Fz)=Jacobi(CG,Grid.Faces[iF],zI,Topo,Grid.Topography,zs[:,:,iF]);
      #(X_Fz1,J_Fz1,dXdx_F1,dXdxI_Fz1,lat_Fz1)=JacobiDG3Neu(CG,Grid.Faces[iF],zI,Topo,Grid.Topography,zs[:,:,iF]);
      X[:,:,:,:,iz,iF]=X_Fz;
      J[:,:,:,iz,iF]=J_Fz;
      dXdx[:,:,:,iz,:,:,iF]=reshape(dXdx_Fz,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
      dXdxI[:,:,:,iz,:,:,iF]=reshape(dXdxI_Fz,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,1,3,3);
      lat[:,:,iF]=lat_Fz;
    end
#   Surface normal
    @views @. FS[:,:,iF] = sqrt(dXdxI[:,:,1,1,3,1,iF] * dXdxI[:,:,1,1,3,1,iF] +
      dXdxI[:,:,1,1,3,2,iF] * dXdxI[:,:,1,1,3,2,iF] + dXdxI[:,:,1,1,3,3,iF] * dXdxI[:,:,1,1,3,3,iF])
    @views @. nS[:,:,1,iF] = dXdxI[:,:,1,1,3,1,iF] / FS[:,:,iF]
    @views @. nS[:,:,2,iF] = dXdxI[:,:,1,1,3,2,iF] / FS[:,:,iF]
    @views @. nS[:,:,3,iF] = dXdxI[:,:,1,1,3,3,iF] / FS[:,:,iF]
  end

  (CG.M,CG.MW,CG.MMass)=MassCG(CG,Global);
  Global.latN=zeros(CG.NumG);
  lat = Global.Metric.lat
  latN = Global.latN
  OP=CG.OrdPoly+1;
  for iF=1:NF
    for jP=1:OP
      for iP=1:OP
        ind=CG.Glob[iP,jP,iF]
        latN[ind] = latN[ind] + lat[iP,jP,iF] * (J[iP,jP,1,1,iF] + J[iP,jP,2,1,iF]) / CG.M[1,ind]
        @views @. dz[:,ind] +=  2.0 * (J[iP,jP,1,1,iF] + J[iP,jP,2,1,iF])^2 / 
          (dXdxI[iP,jP,1,:,3,3,iF] + dXdxI[iP,jP,2,:,3,3,iF])  / CG.M[:,ind]
        @inbounds for iz=1:nz
          if Global.Grid.Form == "Sphere"
            r = norm(0.5 .* (X[iP,jP,1,:,iz,iF] .+ X[iP,jP,2,:,iz,iF]))
            zP[iz,ind] += max(r-Global.Grid.Rad, 0.0) * 
              (J[iP,jP,1,1,iF] + J[iP,jP,2,1,iF]) / CG.M[iz,ind]
          else
            zP[iz,ind] += 0.5 * (X[iP,jP,1,3,iz,iF] + X[iP,jP,2,3,iz,iF]) * 
              (J[iP,jP,1,1,iF] + J[iP,jP,2,1,iF]) / CG.M[iz,ind]
          end
        end
      end
    end
  end
  ExchangeData!(dz,Global.Exchange)
  ExchangeData!(latN,Global.Exchange)
  ExchangeData!(zP,Global.Exchange)
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

