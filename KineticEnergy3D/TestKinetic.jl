using LinearAlgebra
using CairoMakie
using GeometryBasics
using GLMakie
using WriteVTK
include("Grid.jl")
include("GaussLobattoQuad.jl")
include("GaussLegendreQuad.jl")
include("Lagrange.jl")
include("Oro.jl")
include("Interpolation.jl")
include("Jacobi.jl")
include("DSS.jl")
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
include("Cache.jl")
include("../src/Model/Parameters.jl")

function TestKinetic()
  PhysParam = PhysParameters()
  Param = Parameters("HillAgnesiCart")
  Nz = 80
  Nx = 100
  Ny = 2
  OrdPolyX=2
  OrdPolyY=2
  OrdPolyZ=1
  H = 15600.0
  Lx= 40000.0
  Ly= 800.0
  x0 = -Lx/2.0
  y0 = -Ly/2.0

  CacheFcn = Cache(Nx,Ny,Nz,OrdPolyX,OrdPolyY,OrdPolyZ)

  Fe = FeElem(OrdPolyX,OrdPolyY,OrdPolyZ)

  Grid3D = Grid(Nx,Ny,Nz,x0,y0,Lx,Ly,H,Oro,Fe,Param)

  Metric3D = Metric(Nx,Ny,Nz,Grid3D.xP,Grid3D.yP,Grid3D.zP,Fe)

  vtkGrid3D = vtkStruct(Fe,Grid3D.xP,Grid3D.yP,Grid3D.zP)

  RhoPos = 1
  uPos = 2
  vPos = 3
  ThPos = 4
  wPos = 1
  FC = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ,4)
  FF = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,1)
  UC = rand(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ,4)
  UF = rand(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,1)

  @views @. UC[:,:,:,:,:,:,vPos] = 0.0
  @views @. UC[:,:,:,:,:,:,ThPos] = 0.0
  for iy = 1 : Ny
    for j = 1 : OrdPolyY +1
      @views @. UC[:,iy,:,:,j,:,uPos] = UC[:,1,:,:,1,:,uPos]  
      @views @. UC[:,iy,:,:,j,:,RhoPos] = UC[:,1,:,:,1,:,RhoPos]  
      @views @. UF[:,iy,:,:,j,:,wPos] = UF[:,1,:,:,1,:,wPos]  
    end  
  end   
  AverageC!(UC[:,:,:,:,:,:,uPos])
  AverageC!(UC[:,:,:,:,:,:,RhoPos])
  @views wF = UF[:,:,:,:,:,:,wPos]
  @views @. wF[:,:,Nz,:,:,OrdPolyZ+1] = 0.0
  AverageF!(wF)

# with Interpolation for higher order
  @views wF = UF[:,:,:,:,:,:,wPos]
# Lower boundary condition
  uF = zeros(OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1)
  vF = zeros(OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1)
  for ix = 1 : Nx
    for iy = 1 : Ny
      for i = 1 : OrdPolyX + 1  
        for j = 1 : OrdPolyY + 1  
          @views uF[i,j,:] = Fe.IntZC2F * UC[ix,iy,1,i,j,:,uPos]
          @views vF[i,j,:] = Fe.IntZC2F * UC[ix,iy,1,i,j,:,vPos]
        end
      end  
      @views @. wF[ix,iy,1,:,:,1] = -(Metric3D.dXdxI[ix,iy,1,:,:,1,3,1]*uF[:,:,1] +
        Metric3D.dXdxI[ix,iy,1,:,:,1,3,2]*vF[:,:,1]) /
        Metric3D.dXdxI[ix,iy,1,:,:,1,3,3]  
    end    
  end    
  #Average wF
  JXY = zeros(Nx,Ny,OrdPolyX+1,OrdPolyY+1)
  @. JXY = Metric3D.J[:,:,1,:,:,1]
  @views AverageFXY!(wF[:,:,1,:,:,1],JXY)

  Curl!(FC,FF,UC,UF,Metric3D,Fe,PhysParam,CacheFcn)

  @views u = UC[:,:,:,:,:,:,uPos]
  @views Rho = UC[:,:,:,:,:,:,RhoPos]
  @views Fu = FC[:,:,:,:,:,:,uPos]
  @. Fu = u * Fu
  IFu = IntC(Fu,Rho,Metric3D.JC,Fe)


  @show "Ende "

end
