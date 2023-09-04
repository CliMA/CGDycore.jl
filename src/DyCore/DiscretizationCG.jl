mutable struct CGStruct{FT<:Real,
                        AT2<:AbstractArray,
                        IT1<:AbstractArray,
                        IT2<:AbstractArray,
                        IT3<:AbstractArray}
    OrdPoly::Int
    OrdPolyZ::Int
    Glob::IT3
    Stencil::IT2
    NumG::Int
    NumI::Int
    w::Array{FT, 1}
    xw::Array{FT, 1}
    xe::Array{FT, 1}
    IntXE2F::Array{FT, 2}
    xwZ::Array{FT, 1}
    IntZE2F::Array{FT, 2}
    DW::AT2
    DWT::Array{FT, 2}
    DS::AT2
    DST::Array{FT, 2}
    DSZ::Array{FT, 2}
    S::Array{FT, 2}
    M::AT2
    MMass::AT2
    MW::AT2
    Boundary::Array{Int, 1}
    MasterSlave::IT1
end
function CGStruct{FT}(backend) where FT<:Real
  OrdPoly = 0
  OrdPolyZ = 0
  Glob = KernelAbstractions.zeros(backend,Int,0,0,0)
  Stencil = KernelAbstractions.zeros(backend,Int,0,0)
  NumG = 0
  NumI = 0
  w = zeros(FT,0)
  xw = zeros(FT,0)
  xe = zeros(FT,0)
  IntXE2F = zeros(FT,0,0)
  xwZ = zeros(FT,0)
  IntZE2F = zeros(FT,0,0)
  DW = KernelAbstractions.zeros(backend,FT,0,0)
  DWT = zeros(FT,0,0)
  DS = KernelAbstractions.zeros(backend,FT,0,0)
  DST = zeros(FT,0,0)
  DSZ = zeros(FT,0,0)
  S = zeros(FT,0,0)
  M = KernelAbstractions.zeros(backend,FT,0,0)
  MMass = KernelAbstractions.zeros(backend,FT,0,0)
  MW = KernelAbstractions.zeros(backend,FT,0,0)
  Boundary = zeros(Int,0)
  MasterSlave = KernelAbstractions.zeros(backend,Int,0)
 return CGStruct{FT,
                 typeof(DW),
                 typeof(MasterSlave),
                 typeof(Stencil),
                 typeof(Glob)}( 
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


function DiscretizationCG(backend,FT,OrdPoly,OrdPolyZ,Jacobi,Global,zs) 
# Discretization
  Grid = Global.Grid
  OP=OrdPoly+1
  OPZ=OrdPolyZ+1
  nz=Grid.nz
  NF=Grid.NumFaces

  CG = CGStruct{FT}(backend)
  Metric = MetricStruct{FT}(backend,OP,OPZ,Global.Grid.NumFaces,nz)
  CG.OrdPoly=OrdPoly
  CG.OrdPolyZ=OrdPolyZ

  (CG.w,CG.xw)=GaussLobattoQuad(CG.OrdPoly)
  (wZ,CG.xwZ)=GaussLobattoQuad(CG.OrdPolyZ)
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

  (DW,DS)=DerivativeMatrixSingle(CG.OrdPoly)
  CG.DS = KernelAbstractions.zeros(backend,FT,size(DS))
  copyto!(CG.DS,DS)
  CG.DW = KernelAbstractions.zeros(backend,FT,size(DW))
  copyto!(CG.DW,DW)
  CG.DST=DS'
  CG.DWT=DW'

  Q = diagm(CG.w) * DS
  CG.S = Q - Q'
  (DWZ,CG.DSZ)=DerivativeMatrixSingle(CG.OrdPolyZ)
  (Glob,CG.NumG,CG.NumI,Stencil,MasterSlave) =
    NumberingFemCG(Grid,OrdPoly)  


  dXdxI = zeros(3,3,OPZ,OP,OP,nz,NF)
  nS = zeros(OP,OP,3,NF)
  FS = zeros(OP,OP,NF)
  dz = zeros(nz,CG.NumG)
  zP = zeros(nz,CG.NumG)
  J = zeros(OP,OP,OPZ,nz,NF)
  lat = zeros(OP,OP,NF)
  X = zeros(OP,OP,OPZ,3,nz,NF)

  for iF = 1 : NF
    for iz = 1 : nz
      zI = [Grid.z[iz],Grid.z[iz+1]]
      @views (X_Fz,J_Fz,dXdx_Fz,dXdxI_Fz,lat_Fz) = Jacobi(CG,Grid.Faces[iF],zI,Topo,Grid.Topography,zs[:,:,iF])
      @views @. X[:,:,:,:,iz,iF] = X_Fz
      @views @. J[:,:,:,iz,iF] = J_Fz
      @views @. dXdxI[:,:,:,:,:,iz,iF] = dXdxI_Fz
      @views @. lat[:,:,iF] = lat_Fz
      if iz == 1
        #   Surface normal
        @views @. FS[:,:,iF] = sqrt(dXdxI_Fz[3,1,1,:,:] * dXdxI_Fz[3,1,1,:,:] +
          dXdxI_Fz[3,2,1,:,:] * dXdxI_Fz[3,2,1,:,:] + dXdxI_Fz[3,3,1,:,:] * dXdxI_Fz[3,3,1,:,:])
        @views @. nS[:,:,1,iF] = dXdxI_Fz[3,1,1,:,:] / FS[:,:,iF]
        @views @. nS[:,:,2,iF] = dXdxI_Fz[3,2,1,:,:] / FS[:,:,iF]
        @views @. nS[:,:,3,iF] = dXdxI_Fz[3,3,1,:,:] / FS[:,:,iF]
      end
    end
  end

  (M,MW,MMass) = MassCG(CG,J,Glob,Global)
  Global.latN = zeros(CG.NumG)
  latN = Global.latN
  @inbounds for iF = 1 : NF
    @inbounds for jP = 1 : OP
      @inbounds for iP = 1 : OP
        ind = Glob[iP,jP,iF]
        latN[ind] = latN[ind] + lat[iP,jP,iF] * (J[iP,jP,1,1,iF] + J[iP,jP,2,1,iF]) / M[1,ind]
        @inbounds for iz=1:nz
          dz[iz,ind] +=  2.0*(J[iP,jP,1,iz,iF] + J[iP,jP,2,iz,iF])^2 / 
          (dXdxI[3,3,1,iP,jP,iz,iF] + dXdxI[3,3,2,iP,jP,iz,iF])  / M[iz,ind]
        end
        @inbounds for iz=1:nz
          if Global.Grid.Form == "Sphere"
            r = norm(0.5 .* (X[iP,jP,1,:,iz,iF] .+ X[iP,jP,2,:,iz,iF]))
            zP[iz,ind] += max(r-Global.Grid.Rad, 0.0) * 
              (J[iP,jP,1,1,iF] + J[iP,jP,2,1,iF]) / M[iz,ind]
          else
            zP[iz,ind] += 0.5 * (X[iP,jP,1,3,iz,iF] + X[iP,jP,2,3,iz,iF]) * 
              (J[iP,jP,1,1,iF] + J[iP,jP,2,1,iF]) / M[iz,ind]
          end
        end
      end
    end
  end
  ExchangeData!(dz,Global.Exchange)
  ExchangeData!(latN,Global.Exchange)
  ExchangeData!(zP,Global.Exchange)

  copyto!(Metric.dXdxI,dXdxI)
  copyto!(Metric.nS,nS)
  copyto!(Metric.FS,FS)
  Metric.dz = KernelAbstractions.zeros(backend,FT,size(dz))
  copyto!(Metric.dz,dz)
  Metric.zP = KernelAbstractions.zeros(backend,FT,size(zP))
  copyto!(Metric.zP,zP)
  copyto!(Metric.J,J)
  copyto!(Metric.lat,lat)
  copyto!(Metric.X,X)
  CG.M = KernelAbstractions.zeros(backend,FT,size(M))
  copyto!(CG.M,M)
  CG.MMass = KernelAbstractions.zeros(backend,FT,size(MMass))
  copyto!(CG.MMass,MMass)
  CG.MW = KernelAbstractions.zeros(backend,FT,size(MW))
  copyto!(CG.MW,MW)
  CG.Glob = KernelAbstractions.zeros(backend,Int,size(Glob))
  copyto!(CG.Glob,Glob)
  CG.Stencil = KernelAbstractions.zeros(backend,Int,size(Stencil))
  copyto!(CG.Stencil,Stencil)
  CG.MasterSlave = KernelAbstractions.zeros(backend,Int,size(MasterSlave))
  copyto!(CG.MasterSlave,MasterSlave)


# Boundary nodes
  for iF = 1 : NF
    Side = 0
    for iE in Grid.Faces[iF].E
       Side += 1 
       if Grid.Edges[iE].F[1] == 0 || Grid.Edges[iE].F[2] == 0
         if Side == 1
           for i in Glob[1:OP-1,1,iF]   
             push!(CG.Boundary,i)
           end
         elseif Side == 2  
           for i in Glob[OP,1:OP-1,iF]
             push!(CG.Boundary,i)
           end
         elseif Side == 3  
           for i in Glob[2:OP,OP,iF]
             push!(CG.Boundary,i)
           end
         elseif Side == 4  
           for i in Glob[1,2:OP,iF]
             push!(CG.Boundary,i)
           end
         end  
       end  
    end
  end  
  return (CG,Metric,Global)
end

function DiscretizationCG(backend,FT,OrdPoly,OrdPolyZ,Jacobi,Global)
  DiscretizationCG(backend,FT,OrdPoly,OrdPolyZ,Jacobi,Global,
  zeros(OrdPoly+1,OrdPoly+1,Global.Grid.NumFaces))
end  
