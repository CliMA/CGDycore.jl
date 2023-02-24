using LinearAlgebra
using CairoMakie
include("GaussLobattoQuad.jl")
include("GaussLegendreQuad.jl")
include("Lagrange.jl")
include("DLagrange.jl")
include("Oro.jl")
include("JacobiDG2.jl")
include("DSS.jl")
include("DSSF.jl")
include("Curl.jl")
include("Div.jl")
include("Grad.jl")
include("FeElem.jl")
include("IntLGL.jl")
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

function TestModel()
  PhysParam = PhysParameters()
  Param = Parameters("WarmBubble2DXCart")
  Nz = 40
  Nx = 40
  OrdPolyX=3
  OrdPolyZ=2
  H = 10000.0
  hHill = 200
  L= 20000.0

  CacheFcn = Cache(Nx,Nz,OrdPolyX,OrdPolyZ)

  Fe = FeElem(OrdPolyX,OrdPolyZ)

  #Grid 
  xP=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  @views @. xP[1,:,1,:]=0
  dx=L/Nx
  for ix = 1 : Nx
    for iz = 1 : Nz  
      for i = 2 : OrdPolyX + 1
        for j = 1 : OrdPolyZ + 1
          xP[ix,iz,i,j] = xP[ix,iz,1,j] + (1 + Fe.xw[i]) / 2 * dx
        end
      end
    end
    if ix < Nx
      @views @. xP[ix+1,:,1,:] = xP[ix,:,OrdPolyX+1,:]
    end
  end  

  @views @. xP[Nx,:,OrdPolyX+1,:] = L

  zP=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  for ix = 1 : Nx
    for i = 1 : OrdPolyX + 1  
      zP[ix,1,i,1]=Oro(xP[ix,1,i,1],Param)
      dzLoc=(H - zP[ix,1,i,1]) / Nz
      for iz = 1 : Nz
        for j = 1 : OrdPolyZ + 1  
          zP[ix,iz,i,j]=zP[ix,iz,i,1]+(1+Fe.zw[j])/2*dzLoc   
        end  
        if iz < Nz
          zP[ix,iz+1,i,1] = zP[ix,iz,i,OrdPolyZ+1]  
        end  
      end  
      zP[ix,Nz,i,OrdPolyZ+1] = H
    end
  end

  Metric2D = Metric(Nx,Nz,xP,zP,Fe)

  RhoPos = 1
  uPos = 2
  ThPos = 3
  wPos = 1
  UC = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ,3)
  FC = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ,3)
  UF = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1,1)
  FF = zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1,1)
  for ix = 1 : Nx
    for iz = 1 : Nz
      @views fProjectLG!(UC[ix,iz,:,:,uPos],fuVel,Metric2D.X[ix,iz,:,:,:],
        Fe,PhysParam,Param)
      @views fProjectLG!(UC[ix,iz,:,:,RhoPos],fRho,Metric2D.X[ix,iz,:,:,:],
        Fe,PhysParam,Param)
      @views fProjectLG!(UC[ix,iz,:,:,ThPos],fTheta,Metric2D.X[ix,iz,:,:,:],
        Fe,PhysParam,Param)
      @views fProject!(UF[ix,iz,:,:,wPos],fwVel,Metric2D.X[ix,iz,:,:,:],
        PhysParam,Param)
    end  
  end
  @views @. UC[:,:,:,:,ThPos] = UC[:,:,:,:,RhoPos] * UC[:,:,:,:,ThPos]
  @views Plot2DLG(UC[:,:,:,:,ThPos]./UC[:,:,:,:,RhoPos],Fe,"Theta")
  @views Plot2DLG(UC[:,:,:,:,RhoPos],Fe,"Rho")
  @views Plot2DF(UF[:,:,:,:,wPos],Fe,"wEnd")

# with Interpolation for higher order
  @views wF = FF[:,:,:,:,wPos]
# Lower boubdary condition
  uF = zeros(OrdPolyX+1,OrdPolyZ+1)
  for ix = 1 : Nx
    for i = OrdPolyX + 1  
      @views uF[i,:] = Fe.IntZLG2LGL * UC[ix,1,i,:,uPos]
    end  
    @views @. wF[ix,1,:,1] = -Metric2D.dXdxI[ix,1,:,1,2,1]*uF[:,1] /
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

  dtau = 0.1
  IterEnd = 1
  UC_n = similar(UC)
  UF_n = similar(UF)
  for Iter = 1 : IterEnd
    @. UC_n = UC
    @. UF_n = UF
    Fcn!(FC,FF,UC,UF,Metric2D,Fe,PhysParam,CacheFcn)
    @. UC = UC_n + 1.0/3.0 * dtau *FC
    @. UF = UF_n + 1.0/3.0 * dtau *FF
    Fcn!(FC,FF,UC,UF,Metric2D,Fe,PhysParam,CacheFcn)
    @. UC = UC_n + 1.0/2.0 * dtau *FC
    @. UF = UF_n + 1.0/2.0 * dtau *FF
    Fcn!(FC,FF,UC,UF,Metric2D,Fe,PhysParam,CacheFcn)
    @. UC = UC_n + dtau *FC
    @. UF = UF_n + dtau * FF
 #  if mod(Iter,100)
 #    @views Plot2DLG(UC[:,:,:,:,ThPos]./UC[:,:,:,:,RhoPos],Fe,"ThetaEnd")
 #    @views Plot2DLG(UC[:,:,:,:,uPos],Fe,"uEnd")
 #    @views Plot2DF(UF[:,:,:,:,wPos],Fe,"wEnd")
 #  end  
  end  
  @views Plot2DLG(UC[:,:,:,:,ThPos]./UC[:,:,:,:,RhoPos],Fe,"ThetaEnd")
  @views Plot2DLG(UC[:,:,:,:,uPos],Fe,"uEnd")
  @views Plot2DF(UF[:,:,:,:,wPos],Fe,"wEnd")

end
