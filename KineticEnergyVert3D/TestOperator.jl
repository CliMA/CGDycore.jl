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
include("Fun.jl")
include("Pressure.jl")
include("Fcn.jl")
include("PhysParameters.jl")
include("Visualization.jl")
include("Cache.jl")
include("../src/Model/Parameters.jl")

function TestOperator()
  PhysParam = PhysParameters()
  Param = Parameters("WarmBubble2DXCart")
  Nz = 200
  Nx = 200
  Ny = 2
  OrdPolyX=4
  OrdPolyY=1
  OrdPolyZ=2
  H = 1.0
  Lx= 1.0
  Ly= 1.0
  x0 = 0.0
  y0 = 0.0

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
  StepTh = 1
  StepRho = 1
  Stepu = 1
  Stepw = 1
  UC = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ,4)
  FC = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ,4)
  UF = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,1)
  FF = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,1)
  uAdv = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ)
  wAdv = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1)
  Div = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ)
  GradX = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ)
  GradZ = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1)
  for ix = 1 : Nx
    for iy = 1 : Ny
      for iz = 1 : Nz
        @views FunProjectC!(UC[ix,iy,iz,:,:,:,RhoPos],RhoFun,Metric3D.X[ix,iy,iz,:,:,:,:],Fe)
        @views FunProjectC!(UC[ix,iy,iz,:,:,:,uPos],uFun,Metric3D.X[ix,iy,iz,:,:,:,:],Fe)
        @views FunProjectF!(UF[ix,iy,iz,:,:,:,wPos],wFun,Metric3D.X[ix,iy,iz,:,:,:,:])
        @views FunProjectC!(Div[ix,iy,iz,:,:,:],DivFun,Metric3D.X[ix,iy,iz,:,:,:,:],Fe)
        @views FunProjectC!(GradX[ix,iy,iz,:,:,:],GradXKin,Metric3D.X[ix,iy,iz,:,:,:,:],Fe)
        @views FunProjectF!(GradZ[ix,iy,iz,:,:,:],GradZKin,Metric3D.X[ix,iy,iz,:,:,:,:])
        @views FunProjectC!(uAdv[ix,iy,iz,:,:,:],AdvuMom,Metric3D.X[ix,iy,iz,:,:,:,:],Fe)
        @views FunProjectF!(wAdv[ix,iy,iz,:,:,:],AdvwMom,Metric3D.X[ix,iy,iz,:,:,:,:])
      end  
    end  
  end

  Curl!(FC,FF,UC,UF,Metric3D,Fe,PhysParam,CacheFcn)
  for i = 1 : OrdPolyX + 1
    for k = 1 : OrdPolyZ
       @show i,k,FC[20,1,20,i,1,k,uPos],uAdv[20,1,20,i,1,k]
    end
  end
  for i = 1 : OrdPolyX + 1
    for k = 1 : OrdPolyZ + 1
       @show i,k,FF[20,1,20,i,1,k,wPos],wAdv[20,1,20,i,1,k]
    end
  end

  Div!(FC,FF,UC,UF,Metric3D,Fe,PhysParam,CacheFcn)
  for i = 1 : OrdPolyX + 1
    for k = 1 : OrdPolyZ
       @show i,k,FC[20,1,20,i,1,k,RhoPos],Div[20,1,20,i,1,k]
    end
  end

  Grad!(FC,FF,UC,UF,Metric3D,Fe,PhysParam,CacheFcn)
  for i = 1 : OrdPolyX + 1
    for k = 1 : OrdPolyZ
       @show i,k,FC[20,1,20,i,1,k,uPos],GradX[20,1,20,i,1,k]
    end
  end
  for i = 1 : OrdPolyX + 1
    for k = 1 : OrdPolyZ + 1
       @show i,k,FF[20,1,20,i,1,k,wPos],GradZ[20,1,20,i,1,k]
    end
  end
end
