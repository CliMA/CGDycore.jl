using LinearAlgebra
using CairoMakie
using GeometryBasics
using GLMakie
include("GaussLobattoQuad.jl")
include("GaussLegendreQuad.jl")
include("Lagrange.jl")
include("Oro.jl")
include("JacobiDG2.jl")
include("DSS.jl")
include("DSSF.jl")
include("Curl.jl")
include("Div.jl")
include("Grad.jl")
include("FeElem.jl")
include("Int.jl")
include("Metric.jl")
include("ModelFun.jl")
include("Pressure.jl")
include("Buoyancy.jl")
include("Fcn.jl")
include("PhysParameters.jl")
include("Visualization.jl")
include("../src/Model/Parameters.jl")

mutable struct Cache
  pC::Array{Float64, 4}
  pF::Array{Float64, 4}
  uF::Array{Float64, 4}
  RhoF::Array{Float64, 4}
  RhoThetaF::Array{Float64, 4}
  FuF::Array{Float64, 4}
  FRhoF::Array{Float64, 4}
  FRhoThetaF::Array{Float64, 4}
  KinF::Array{Float64, 4}
  KinC::Array{Float64, 4}
  Column::Array{Float64, 3}
  Block::Array{Float64, 3}
end

function Cache(Nx,Nz,OrdPolyX,OrdPolyZ)
  pC = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ)
  pF = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  uF = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  RhoF = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  RhoThetaF = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  FuF = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  FRhoF = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  FRhoThetaF = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  KinF = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  KinC = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ)
  Column = zeros(OrdPolyX+1,OrdPolyZ+1,Nz)
  Block = zeros(OrdPolyX+1,OrdPolyZ+1,6)

  return Cache(
    pC,
    pF,
    uF,
    RhoF,
    RhoThetaF,
    FuF,
    FRhoF,
    FRhoThetaF,
    KinF,
    KinC,
    Column,
    Block,
  )
end  

function TestKinetic()

PhysParam = PhysParameters()
  Param = Parameters("HillAgnesiCart")
  Nz = 40
  Nx = 100
  OrdPolyX=4
  OrdPolyZ=1
  H = 15600.0
  Lx= 60000.0
  x0 = -Lx/2.0

  CacheFcn = Cache(Nx,Nz,OrdPolyX,OrdPolyZ)

  Fe = FeElem(OrdPolyX,OrdPolyZ)

  Grid2D = Grid(Nx,Nz,x0,Lx,H,Oro,Fe,Param)

  Metric2D = Metric(Nx,Nz,Grid2D.xP,Grid2D.zP,Fe)

  RhoPos = 1 
  uPos = 2 
  ThPos = 3 
  wPos = 1 

  u=rand(Nx,Nz,OrdPolyX+1,OrdPolyZ)
  wF=2.0*rand(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  Rho=rand(Nx,Nz,OrdPolyX+1,OrdPolyZ)
  @. Rho = Rho + 1.0
  Average!(u)
  Average!(Rho)
  @views @. wF[:,Nz,:,OrdPolyZ+1] = 0.0
  AverageF!(wF)

# with Interpolation for higher order
  Fu = similar(u)
  FRho = similar(u)
  FwF = similar(wF)
  FuF = similar(wF)
  FRhoF = similar(wF)
  uF = similar(wF)
  RhoF = similar(wF)
  KinF = similar(wF)
  Kin = similar(u)
  for ix = 1 : Nx
    for iz = 1 : Nz
      for i = 1 : OrdPolyX +1
        @views uF[ix,iz,i,:] = Fe.IntZC2F * u[ix,iz,i,:]
        @views RhoF[ix,iz,i,:] = Fe.IntZC2F * Rho[ix,iz,i,:]
      end
    end
  end
  for ix = 1 : Nx
    @views @. wF[ix,1,:,1] = -Metric2D.dXdxI[ix,1,:,1,2,1]*uF[ix,1,:,1] /
      Metric2D.dXdxI[ix,1,:,1,2,2]
  end    
  #Average wF
  for ix = 2 : Nx
    ww = 0.5 * ( wF[ix-1,1,OrdPolyX+1,1] + wF[ix,1,1,1]) 
    wF[ix-1,1,OrdPolyX+1,1] = ww
    wF[ix,1,1,1] = ww
  end
  ww = 0.5 * ( wF[Nx,1,OrdPolyX+1,1] + wF[1,1,1,1])
  wF[Nx,1,OrdPolyX+1,1] = ww
  wF[1,1,1,1] = ww
  @. KinF = 0.5 * (wF * wF + uF * uF)
  for ix = 1 : Nx
    for iz = 1 : Nz
      for i = 1 : OrdPolyX +1
        @views Kin[ix,iz,i,:] = Fe.IntZF2C * KinF[ix,iz,i,:]
        @views KinF[ix,iz,i,:] = Fe.IntZC2F * Kin[ix,iz,i,:]
      end
    end
  end

  @. FuF = 0.0
  @. FwF = 0.0
  @. FRhoF = 0.0
  Curl!(FuF,FwF,uF,wF,RhoF,Fe,Metric2D,CacheFcn)
  Div!(FRhoF,uF,wF,RhoF,Fe,Metric2D)
  RhoGrad!(FuF,FwF,KinF,RhoF,Fe,Metric2D)
  su = IntF(FuF.*uF,Fe)
  sw = IntF(FwF.*wF,Fe)
  @show su
  @show sw
  sRho = IntF(FRhoF.*KinF,Fe)
  @show sRho


  DSSF!(FwF,RhoF,Metric2D.J)
  DSS!(FuF,RhoF,Metric2D.J)
  DSS!(FRhoF,Metric2D.J)
  for ix = 1 : Nx
    for iz = 1 : Nz
      for i = 1 : OrdPolyX +1
        @views Fu[ix,iz,i,:] =  Fe.P * FuF[ix,iz,i,:]
        @views FRho[ix,iz,i,:] = Fe.P * FRhoF[ix,iz,i,:]
      end
    end
  end
  @. FwF = wF * FwF
  IFwF = IntF(FwF,RhoF,Metric2D.J,Fe)
  @. Fu = u * Fu
  IFu = IntC(Fu,Rho,Metric2D.JC,Fe)
  @. FRho = FRho * Kin
  IFRho = IntC(FRho,Metric2D.JC,Fe)
  @show IFu , IFwF , IFRho
  @show IFu + IFwF + IFRho

end
