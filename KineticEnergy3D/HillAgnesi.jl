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
include("Cache.jl")
include("../src/Model/Parameters.jl")

function HillAgnesi()
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
  UC = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ,4)
  FC = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ,4)
  UF = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,1)
  FF = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,1)
  for ix = 1 : Nx
    for iy = 1 : Ny
      for iz = 1 : Nz
        @views fProjectC!(UC[ix,iy,iz,:,:,:,uPos],fuVel,Metric3D.X[ix,iy,iz,:,:,:,:],
          Fe,PhysParam,Param)
        @views fProjectC!(UC[ix,iy,iz,:,:,:,vPos],fvVel,Metric3D.X[ix,iy,iz,:,:,:,:],
          Fe,PhysParam,Param)
        @views fProjectC!(UC[ix,iy,iz,:,:,:,RhoPos],fRho,Metric3D.X[ix,iy,iz,:,:,:,:],
          Fe,PhysParam,Param)
        @views fProjectC!(UC[ix,iy,iz,:,:,:,ThPos],fTheta,Metric3D.X[ix,iy,iz,:,:,:,:],
          Fe,PhysParam,Param)
        @views fProjectF!(UF[ix,iy,iz,:,:,:,wPos],fwVel,Metric3D.X[ix,iy,iz,:,:,:,:],
          PhysParam,Param)
      end  
    end  
  end
  @views @. UC[:,:,:,:,:,:,ThPos] = UC[:,:,:,:,:,:,RhoPos] * UC[:,:,:,:,:,:,ThPos]

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

  @views vtkPlot3DC(UC[:,:,:,:,:,:,RhoPos],Fe,vtkGrid3D,"Rho")
  @views vtkPlot3DC(UC[:,:,:,:,:,:,uPos],Fe,vtkGrid3D,"uVel")
  @views vtkPlot3DF(UF[:,:,:,:,:,:,wPos],Fe,vtkGrid3D,"wVel")
  @views vtkPlot3DC(UC[:,:,:,:,:,:,ThPos]./UC[:,:,:,:,:,:,RhoPos],Fe,vtkGrid3D,"Theta")
  vtkGrid3D.Step +=1

  dtau = 0.1
  IterEnd = 4000 #10000
  UC_n = similar(UC)
  UF_n = similar(UF)
  for Iter = 1 : IterEnd
    @show Iter  
    @. UC_n = UC
    @. UF_n = UF
    Fcn!(FC,FF,UC,UF,Metric3D,Fe,PhysParam,CacheFcn)
    @. UC = UC_n + 1.0/3.0 * dtau *FC
    @. UF = UF_n + 1.0/3.0 * dtau *FF
    Fcn!(FC,FF,UC,UF,Metric3D,Fe,PhysParam,CacheFcn)
    @. UC = UC_n + 1.0/2.0 * dtau *FC
    @. UF = UF_n + 1.0/2.0 * dtau *FF
    Fcn!(FC,FF,UC,UF,Metric3D,Fe,PhysParam,CacheFcn)
    @. UC = UC_n + dtau *FC
    @. UF = UF_n + dtau * FF
    if mod(Iter,400) == 0
      @views vtkPlot3DC(UC[:,:,:,:,:,:,RhoPos],Fe,vtkGrid3D,"Rho")
      @views vtkPlot3DC(UC[:,:,:,:,:,:,uPos],Fe,vtkGrid3D,"uVel")
      @views vtkPlot3DF(UF[:,:,:,:,:,:,wPos],Fe,vtkGrid3D,"wVel")
      @views vtkPlot3DC(UC[:,:,:,:,:,:,ThPos]./UC[:,:,:,:,:,:,RhoPos],Fe,vtkGrid3D,"Theta")
      vtkGrid3D.Step +=1
    end  
  end  
  @views vtkPlot3DC(UC[:,:,:,:,:,:,RhoPos],Fe,vtkGrid3D,"Rho")
  @views vtkPlot3DC(UC[:,:,:,:,:,:,uPos],Fe,vtkGrid3D,"uVel")
  @views vtkPlot3DF(UF[:,:,:,:,:,:,wPos],Fe,vtkGrid3D,"wVel")
  @views vtkPlot3DC(UC[:,:,:,:,:,:,ThPos]./UC[:,:,:,:,:,:,RhoPos],Fe,vtkGrid3D,"Theta")

  @show "Ende "

end
