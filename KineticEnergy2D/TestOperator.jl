using LinearAlgebra
include("GaussLobattoQuad.jl")
include("GaussLegendreQuad.jl")
include("Lagrange.jl")
include("DLagrange.jl")
include("Oro.jl")
include("JacobiDG2.jl")
include("DSS.jl")
include("DSSF.jl")
include("udwdx.jl")
include("wdwdx.jl")
include("IntCell.jl")
include("IntFace.jl")
include("dcdx.jl")
include("divdx.jl")
include("Sdwdz.jl")
include("Sdudz.jl")
include("dRhoSdz.jl")
include("Curl.jl")
include("Div.jl")
include("Grad.jl")
include("FeElem.jl")
include("IntLGL.jl")
include("Metric.jl")
include("Fun.jl")

function TestOperator()
  Nz = 80
  Nx = 80
  OrdPolyX=4
  OrdPolyZ=4
  H = 1 #400
  hHill = 200
  L= 1 #1000


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
#     zP[ix,1,i,1]=Oro(xP[ix,1,i,1],L,hHill)
      zP[ix,1,i,1] = 0.0
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

  u=rand(Nx,Nz,OrdPolyX+1,OrdPolyZ)
  wF=2.0*rand(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  Rho=rand(Nx,Nz,OrdPolyX+1,OrdPolyZ)
  @. Rho = Rho + 1.0
  RhoTh=rand(Nx,Nz,OrdPolyX+1,OrdPolyZ)
  Div=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ)
  GradZ=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  GradX=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ)
  uAdv=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ)
  wAdv=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
# @. u = 1.0
  for ix = 1 : Nx
    for iz = 1 : Nz
      @views FunProjectLG!(u[ix,iz,:,:],uFun,Metric2D.X[ix,iz,:,:,:],Fe)  
      @views FunProjectLG!(Rho[ix,iz,:,:],RhoFun,Metric2D.X[ix,iz,:,:,:],Fe)  
      @views FunProject!(wF[ix,iz,:,:],wFun,Metric2D.X[ix,iz,:,:,:])  
      @views FunProjectLG!(Div[ix,iz,:,:],DivFun,Metric2D.X[ix,iz,:,:,:],Fe)  
      @views FunProjectLG!(GradX[ix,iz,:,:],GradXKin,Metric2D.X[ix,iz,:,:,:],Fe)  
      @views FunProject!(GradZ[ix,iz,:,:],GradZKin,Metric2D.X[ix,iz,:,:,:])  
      @views FunProjectLG!(uAdv[ix,iz,:,:],AdvuMom,Metric2D.X[ix,iz,:,:,:],Fe)  
      @views FunProject!(wAdv[ix,iz,:,:],AdvwMom,Metric2D.X[ix,iz,:,:,:])  
    end  
  end  

# with Interpolation for higher order
  Fu = similar(u)
  FRho = similar(u)
  FRhoTh = similar(u)
  FwF = similar(wF)
  FuLGL = similar(wF)
  FRhoLGL = similar(wF)
  FRhoThLGL = similar(wF)
  uLGL = similar(wF)
  RhoLGL = similar(wF)
  RhoThLGL = similar(wF)
  KinLGL = similar(wF)
  Kin = similar(u)
  pLGL = similar(wF)
  p = similar(u)
  for ix = 1 : Nx
    for iz = 1 : Nz
      for i = 1 : OrdPolyX +1
        @views uLGL[ix,iz,i,:] = Fe.IntZLG2LGL * u[ix,iz,i,:]
        @views RhoLGL[ix,iz,i,:] = Fe.IntZLG2LGL * Rho[ix,iz,i,:]
        @views RhoThLGL[ix,iz,i,:] = Fe.IntZLG2LGL * RhoTh[ix,iz,i,:]
        @views pLGL[ix,iz,i,:] = Fe.IntZLG2LGL * p[ix,iz,i,:]
      end
    end
  end
  @. KinLGL = 0.5 * (wF * wF + uLGL * uLGL)
  for ix = 1 : Nx
    for iz = 1 : Nz
      for i = 1 : OrdPolyX +1
        @views Kin[ix,iz,i,:] = Fe.IntZLGL2LG * KinLGL[ix,iz,i,:]
        @views KinLGL[ix,iz,i,:] = Fe.IntZLG2LGL * Kin[ix,iz,i,:]
      end
    end
  end

  @. FuLGL = 0.0
  @. FwF = 0.0
  Curl!(FuLGL,FwF,uLGL,wF,RhoLGL,Fe,Metric2D)
  DSSF!(FwF,RhoLGL,Metric2D.J)
  DSS!(FuLGL,RhoLGL,Metric2D.J)
  for ix = 1 : Nx
    for iz = 1 : Nz
      for i = 1 : OrdPolyX +1
        @views Fu[ix,iz,i,:] =  Fe.P * FuLGL[ix,iz,i,:]
      end
    end
  end
  for i = 1 : OrdPolyX + 1
    for j = 1 : OrdPolyZ + 1
       @show i,j,FwF[20,20,i,j],wAdv[20,20,i,j] 
    end   
  end   
  for i = 1 : OrdPolyX + 1
    for j = 1 : OrdPolyZ
       @show i,j,Fu[20,20,i,j],uAdv[20,20,i,j] 
    end   
  end   


  @. FuLGL = 0.0
  @. FwF = 0.0
  @. FRhoLGL = 0.0
  @. FRhoThLGL = 0.0
  Div!(FRhoLGL,uLGL,wF,RhoLGL,Fe,Metric2D)
  RhoGrad!(FuLGL,FwF,KinLGL,RhoLGL,Fe,Metric2D)

  #DivTh!(FRhoThLGL,RhoThLGL,uLGL,wF,RhoLGL,Fe,Metric2D)

  DSSF!(FwF,RhoLGL,Metric2D.J)
  DSS!(FuLGL,RhoLGL,Metric2D.J)
  DSS!(FRhoLGL,Metric2D.J)
  for ix = 1 : Nx
    for iz = 1 : Nz
      for i = 1 : OrdPolyX +1
        @views Fu[ix,iz,i,:] =  Fe.P * FuLGL[ix,iz,i,:]
        @views FRho[ix,iz,i,:] = Fe.P * FRhoLGL[ix,iz,i,:]
      end
    end
  end
  for i = 1 : OrdPolyX + 1
    for j = 1 : OrdPolyZ + 1
       @show i,j,FwF[20,20,i,j],GradZ[20,20,i,j] 
    end   
  end   
  for i = 1 : OrdPolyX + 1
    for j = 1 : OrdPolyZ
       @show i,j,Fu[20,20,i,j],GradX[20,20,i,j] 
    end   
  end   
  for i = 1 : OrdPolyX + 1
    for j = 1 : OrdPolyZ
       @show i,j,FRho[20,20,i,j],Div[20,20,i,j]
    end   
  end   
end
