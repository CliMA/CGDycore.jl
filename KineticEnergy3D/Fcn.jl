function Curl!(FC,FF,UC,UF,Metric,Fe,PhysParam,Cache)
  RhoCPos = 1
  uCPos = 2
  vCPos = 3
  wFPos = 1

  Nx = size(FC,1)
  Ny = size(FC,2)
  Nz = size(FC,3)

  @views RhoC = UC[:,:,:,:,:,1]
  @views uC =   UC[:,:,:,:,:,2]
  @views vC =   UC[:,:,:,:,:,3]
  @views FuC =  FC[:,:,:,:,:,2]
  @views FvC =  FC[:,:,:,:,:,3]
  @views wF =   UF[:,:,:,:,:,wFPos]
  @views FwF =  FF[:,:,:,:,:,wFPos]

  @views JRhoC = Cache.KinC
  @views JRhoF = Cache.KinF


  @. FuC = 0.0
  @. FvC = 0.0
  @. FwF = 0.0
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @views CurlColumn!(FuC[ix,iy,:,:,:],FvC[ix,iy,:,:,:],FwF[ix,iy,:,:,:],
        uC[ix,iy,:,:,:],vC[ix,iy,:,:,:],wF[ix,iy,:,:,:],
        RhoC[ix,iy,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
    end
  end
  @views @. JRhoF[:,:,1,:,:] = Metric.J[:,:,1,:,:,1] * RhoC[:,:,1,:,:]
  @views @. JRhoF[:,:,2:Nz,:,:] = (Metric.J[:,:,1:Nz-1,:,:,2] * RhoC[:,:,1:Nz-1,:,:] +
    Metric.J[:,:,2:Nz,:,:,1] * RhoC[:,:,2:Nz,:,:])                                     
  @views @. JRhoF[:,:,Nz,:,:] = Metric.J[:,:,Nz,:,:,2] * RhoC[:,:,Nz,:,:] 
  DSSF!(FwF,JRhoF)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(FuC,JRhoC)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(FvC,JRhoC)
end

function Div!(FC,FF,UC,UF,Metric,Fe,PhysParam,Cache)
  RhoCPos = 1
  uCPos = 2
  vCPos = 3
  wFPos = 1

  Nx = size(FC,1)
  Ny = size(FC,2)
  Nz = size(FC,3)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  @views RhoC = UC[:,:,:,:,:,1]
  @views uC =   UC[:,:,:,:,:,2]
  @views vC =   UC[:,:,:,:,:,3]
  @views FRhoC =  FC[:,:,:,:,:,1]

  @views wF = UF[:,:,:,:,:,wFPos]

  @views JRhoC = Cache.KinC

  JRhoF = Cache.KinF

  @. FRhoC = 0.0
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @views DivRhoColumn!(FRhoC[ix,iy,:,:,:],uC[ix,iy,:,:,:],vC[ix,iy,:,:,:],
        wF[ix,iy,:,:,:],RhoC[ix,iy,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
    end
  end
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) 
  DSSC!(FRhoC,JRhoC)
end

function Grad!(FC,FF,UC,UF,Metric,Fe,PhysParam,Cache)
  RhoCPos = 1
  uCPos = 2
  vCPos = 3
  wFPos = 1

  Nx = size(FC,1)
  Ny = size(FC,2)
  Nz = size(FC,3)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  @views RhoC = UC[:,:,:,:,:,1]
  @views uC =   UC[:,:,:,:,:,2]
  @views vC =   UC[:,:,:,:,:,3]
  @views FuC =  FC[:,:,:,:,:,2]
  @views FvC =  FC[:,:,:,:,:,3]
  @views wF =   UF[:,:,:,:,:,wFPos]
  @views FwF =  FF[:,:,:,:,:,wFPos]
  KinC = Cache.KinC
  JRhoC = Cache.KinC
  JRhoF = Cache.KinF

  @views @. KinC = uC * uC + vC * vC + 0.5 *(wF[:,:,1:Nz,:,:,] * wF[:,:,1:Nz,:,:] +
    wF[:,:,2:Nz+1,:,:] * wF[:,:,2:Nz+1,:,:])
  @. FuC = 0.0
  @. FvC = 0.0
  @. FwF = 0.0
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @views RhoGradColumn!(FuC[ix,iy,:,:,:],FvC[ix,iy,:,:,:],FwF[ix,iy,:,:,:],
        KinC[ix,iy,:,:,:,:],RhoC[ix,iy,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
    end
  end
  @views @. JRhoF[:,:,1,:,:] = Metric.J[:,:,1,:,:,1] * RhoC[:,:,1,:,:]
  @views @. JRhoF[:,:,2:Nz,:,:] = (Metric.J[:,:,1:Nz-1,:,:,2] * RhoC[:,:,1:Nz-1,:,:] +
    Metric.J[:,:,2:Nz,:,:,1] * RhoC[:,:,2:Nz,:,:])
  @views @. JRhoF[:,:,Nz,:,:] = Metric.J[:,:,Nz,:,:,2] * RhoC[:,:,Nz,:,:]
  DSSF!(FwF,JRhoF)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(FuC,JRhoC)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(FvC,JRhoC)

end

function Hydrostatic!(FC,FF,UC,UF,Metric,Fe,PhysParam,Cache)

  IntZC2F = Fe.IntZC2F
  IntZF2C = Fe.IntZF2C
  P = Fe.P

  RhoCPos = 1
  uCPos = 2
  vCPos = 3
  ThCPos = 4
  wFPos = 1

  Grav = PhysParam.Grav

  Nx = size(FC,1)
  Ny = size(FC,2)
  Nz = size(FC,3)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  pC = Cache.pC
  pF = Cache.pF
  @views RhoF = Cache.UUF[:,:,:,:,:,:,1]
  @views uF = Cache.UUF[:,:,:,:,:,:,2]
  @views vF = Cache.UUF[:,:,:,:,:,:,3]
  @views RhoThetaF = Cache.UUF[:,:,:,:,:,:,4]
  @views FRhoF = Cache.FUUF[:,:,:,:,:,:,1]
  @views FuF = Cache.FUUF[:,:,:,:,:,:,2]
  @views FvF = Cache.FUUF[:,:,:,:,:,:,3]
  @views FRhoThetaF = Cache.FUUF[:,:,:,:,:,:,4]
  @views FRhoThetaF = Cache.FUUF[:,:,:,:,:,:,4]
  JRhoF = Cache.KinF
  PhiF = Cache.KinF

  @views wF = UF[:,:,:,:,:,:,wFPos]
  @views FwF = FF[:,:,:,:,:,:,wFPos]

  BoundaryW!(wF,UC,Metric.dXdxI,Fe,Cache)
  @. wF[:,:,Nz,:,:,OrdPolyZ+1] = 0.0

  @views Pressure!(pC,UC[:,:,:,:,:,:,ThCPos],PhysParam)
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @inbounds for iz = 1 : Nz
        @inbounds for i = 1 : OrdPolyX +1
          @inbounds for j = 1 : OrdPolyY +1
            @inbounds for k = 1 : OrdPolyZ + 1
              uF[ix,iy,iz,i,j,k] = 0.0
              vF[ix,iy,iz,i,j,k] = 0.0
              RhoF[ix,iy,iz,i,j,k] = 0.0
              RhoThetaF[ix,iy,iz,i,j,k] = 0.0
              pF[ix,iy,iz,i,j,k] = 0.0
              @inbounds for l = 1 : OrdPolyZ 
                uF[ix,iy,iz,i,j,k] += IntZC2F[k,l] * UC[ix,iy,iz,i,j,l,uCPos]
                vF[ix,iy,iz,i,j,k] += IntZC2F[k,l] * UC[ix,iy,iz,i,j,l,vCPos]
                RhoF[ix,iy,iz,i,j,k] += IntZC2F[k,l] * UC[ix,iy,iz,i,j,l,RhoCPos]
                RhoThetaF[ix,iy,iz,i,j,k] += IntZC2F[k,l] * UC[ix,iy,iz,i,j,l,ThCPos]
                pF[ix,iy,iz,i,j,k] += IntZC2F[k,l] * pC[ix,iy,iz,i,j,l]
              end
            end  
          end
        end
      end
    end
  end

  @. FuF = 0.0
  @. FvF = 0.0
  @. FwF = 0.0
  @. FRhoF = 0.0
  @. FRhoThetaF = 0.0
  @views @. PhiF = Grav * Metric.X[:,:,:,:,:,:,3]
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
#     @views RhoGradColumn!(FuF[ix,iy,:,:,:,:],FvF[ix,iy,:,:,:,:],FwF[ix,iy,:,:,:,:],
#       PhiF[ix,iy,:,:,:,:],RhoF[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
      @views GradColumn!(FuF[ix,iy,:,:,:,:],FvF[ix,iy,:,:,:,:],FwF[ix,iy,:,:,:,:],
        pF[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
      @views  Buoyancy!(FwF[ix,iy,:,:,:,:],RhoF[ix,iy,:,:,:,:],Metric.J[ix,iy,:,:,:,:],PhysParam)
    end
  end  

  @. JRhoF = Metric.J * RhoF
  DSSF!(FwF,JRhoF)
  @. JRhoF = Metric.J * RhoF
  DSSC!(FuF,JRhoF)
  @. JRhoF = Metric.J * RhoF
  DSSC!(FvF,JRhoF)
  @. JRhoF = Metric.J 
  DSSC!(FRhoF,JRhoF)
  @. JRhoF = Metric.J 
  DSSC!(FRhoThetaF,JRhoF)
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @inbounds for iz = 1 : Nz
        @inbounds for i = 1 : OrdPolyX +1
          @inbounds for j = 1 : OrdPolyY +1
            @inbounds for k = 1 : OrdPolyZ 
              FC[ix,iy,iz,i,j,k,uCPos] = 0.0
              FC[ix,iy,iz,i,j,k,vCPos] = 0.0
              FC[ix,iy,iz,i,j,k,RhoCPos] = 0.0
              FC[ix,iy,iz,i,j,k,ThCPos] = 0.0
              @inbounds for l = 1 : OrdPolyZ +1
                FC[ix,iy,iz,i,j,k,uCPos] +=  P[k,l] * FuF[ix,iy,iz,i,j,l]
                FC[ix,iy,iz,i,j,k,vCPos] +=  P[k,l] * FvF[ix,iy,iz,i,j,l]
                FC[ix,iy,iz,i,j,k,RhoCPos] += P[k,l] * FRhoF[ix,iy,iz,i,j,l]
                FC[ix,iy,iz,i,j,k,ThCPos] +=  P[k,l] * FRhoThetaF[ix,iy,iz,i,j,l]
              end
            end  
          end
        end
      end
    end
  end
  @show size(FwF)
  for k = 1 : OrdPolyZ +1
    @show FwF[10,1,4,2,1,k],pF[10,1,4,2,1,k]
  end
  for k = 1 : OrdPolyZ +1
    @show FwF[10,1,5,2,1,k],pF[10,1,5,2,1,k]
  end
end  

function Fcn!(FC,FF,UC,UF,Metric,Fe,PhysParam,Cache,Koeff)

  IntZC2F = Fe.IntZC2F
  IntZF2C = Fe.IntZF2C
  P = Fe.P

  RhoCPos = 1
  uCPos = 2
  vCPos = 3
  ThCPos = 4
  wFPos = 1

  Nx = size(FC,1)
  Ny = size(FC,2)
  Nz = size(FC,3)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  @views RhoC = UC[:,:,:,:,:,RhoCPos]
  @views uC =   UC[:,:,:,:,:,uCPos]
  @views vC =   UC[:,:,:,:,:,vCPos]
  @views RhoThetaC = UC[:,:,:,:,:,ThCPos]
  @views FRhoC =  FC[:,:,:,:,:,RhoCPos]
  @views FuC =    FC[:,:,:,:,:,uCPos]
  @views FvC =  FC[:,:,:,:,:,vCPos]
  @views FRhoThetaC =  FC[:,:,:,:,:,ThCPos]
  @views wF =   UF[:,:,:,:,:,wFPos]
  @views FwF =  FF[:,:,:,:,:,wFPos]
  KinC = Cache.KinC
  JRhoC = Cache.KinC
  JRhoF = Cache.KinF
  pC = Cache.pC[:,:,:,:,:,1]
  @views CuC = Cache.CF[:,:,:,:,:,1,1]
  @views CvC = Cache.CF[:,:,:,:,:,1,2]
  @views GuC = Cache.CF[:,:,:,:,:,1,3]
  @views GvC = Cache.CF[:,:,:,:,:,1,4]
  @views DThC = Cache.CF[:,:,:,:,:,1,5]

  BoundaryW!(wF,UC,Metric.dXdxI,Fe,Cache)
  @. wF[:,:,Nz+1,:,:] = 0.0

  @views Pressure!(pC,RhoThetaC,PhysParam)

  @. CuC = 0.0
  @. CvC = 0.0
  @. GuC = 0.0
  @. GuC = 0.0
  @. DThC = 0.0
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @views RotCurl!(CuC[ix,iy,:,:,:],CvC[ix,iy,:,:,:],
        uC[ix,iy,:,:,:],vC[ix,iy,:,:,:],Fe,
        Metric.dXdxI[ix,iy,:,:,:,:,:,:],Metric.J[ix,iy,:,:,:,:],Cache)
      @views GradDiv!(GuC[ix,iy,:,:,:,:],GvC[ix,iy,:,:,:,:],
        uC[ix,iy,:,:,:,:],vC[ix,iy,:,:,:,:],Fe,
        Metric.dXdxI[ix,iy,:,:,:,:,:,:],Metric.J[ix,iy,:,:,:,:],Cache)
      @views FDivRhoGrad!(DThC[ix,iy,:,:,:,:],RhoThetaC[ix,iy,:,:,:,:],
        RhoC[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],
        Metric.J[ix,iy,:,:,:,:],Cache)
    end
  end  
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(CuC,JRhoC)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(CvC,JRhoC)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(GuC,JRhoC)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(GvC,JRhoC)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(DThC,JRhoC)

  @views @. KinC = uC * uC + vC * vC + 0.5 *(wF[:,:,1:Nz,:,:,] * wF[:,:,1:Nz,:,:] +
    wF[:,:,2:Nz+1,:,:] * wF[:,:,2:Nz+1,:,:])

  @. FuC = 0.0
  @. FvC = 0.0
  @. FwF = 0.0
  @. FRhoC = 0.0
  @. FRhoThetaC = 0.0
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @views RotCurl!(FuC[ix,iy,:,:,:],FvC[ix,iy,:,:,:],
        CuC[ix,iy,:,:,:],CvC[ix,iy,:,:,:],RhoC[ix,iy,:,:,:],Fe,
        Metric.dXdxI[ix,iy,:,:,:,:,:,:],Metric.J[ix,iy,:,:,:,:],Cache,Koeff)
      @views GradDiv!(FuC[ix,iy,:,:,:],FvC[ix,iy,:,:,:],
        GuC[ix,iy,:,:,:],GvC[ix,iy,:,:,:],RhoC[ix,iy,:,:,:],Fe,
        Metric.dXdxI[ix,iy,:,:,:,:,:,:],Metric.J[ix,iy,:,:,:,:],Cache,Koeff)
      @views FDivRhoGrad!(FRhoThetaC[ix,iy,:,:,:],DThC[ix,iy,:,:,:,:],
        RhoC[ix,iy,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],
        Metric.J[ix,iy,:,:,:,:],Cache,Koeff)
#     @views CurlColumn!(FuC[ix,iy,:,:,:],FvC[ix,iy,:,:,:],FwF[ix,iy,:,:,:],
#       uC[ix,iy,:,:,:],vC[ix,iy,:,:,:],wF[ix,iy,:,:,:],
#       RhoC[ix,iy,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
#     @views RhoGradColumn!(FuC[ix,iy,:,:,:],FvC[ix,iy,:,:,:],FwF[ix,iy,:,:,:],
#       KinC[ix,iy,:,:,:],RhoC[ix,iy,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
      @views GradColumn!(FuC[ix,iy,:,:,:],FvC[ix,iy,:,:,:],FwF[ix,iy,:,:,:],
        pC[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
      @views DivRhoColumn!(FRhoC[ix,iy,:,:,:],uC[ix,iy,:,:,:],vC[ix,iy,:,:,:],
        wF[ix,iy,:,:,:,:],RhoC[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
      @views DivRhoThetaColumn!(FRhoThetaC[ix,iy,:,:,:],uC[ix,iy,:,:,:],vC[ix,iy,:,:,:],
        wF[ix,iy,:,:,:,:],RhoThetaC[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
      @views  Buoyancy!(FwF[ix,iy,:,:,:],RhoC[ix,iy,:,:,:],Metric.J[ix,iy,:,:,:,:],PhysParam)
    end
  end  
  @views @. JRhoF[:,:,1,:,:] = Metric.J[:,:,1,:,:,1] * RhoC[:,:,1,:,:]
  @views @. JRhoF[:,:,2:Nz,:,:] = (Metric.J[:,:,1:Nz-1,:,:,2] * RhoC[:,:,1:Nz-1,:,:] +
    Metric.J[:,:,2:Nz,:,:,1] * RhoC[:,:,2:Nz,:,:])
  @views @. JRhoF[:,:,Nz+1,:,:] = Metric.J[:,:,Nz,:,:,2] * RhoC[:,:,Nz,:,:]
  DSSF!(FwF,JRhoF)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2])
  DSSC!(FRhoC,JRhoC)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(FuC,JRhoC)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(FvC,JRhoC)
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2])
  DSSC!(FRhoThetaC,JRhoC)

  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @views Damping!(FwF[ix,iy,:,:,:],wF[ix,iy,:,:,:],Metric.X[ix,iy,:,:,:,:,:],Fe)
    end
  end  
end  

function Advection!(FC,c,uC,vC,Fe,Metric,Cache)

  Nx = size(FC,1)
  Ny = size(FC,2)
  Nz = size(FC,3)

  @views JRhoC = Cache.KinC

  @views @.  FC = 0.0
  for ix = 1 : Nx
    for iy = 1 : Ny
      for iz = 1 : Nz
        @views AdvectionBlock!(FC[ix,iy,iz,:,:],c[ix,iy,iz,:,:],uC[ix,iy,iz,:,:],vC[ix,iy,iz,:,:],Fe,
          Metric.dXdxI[ix,iy,iz,:,:,:,:,:],Metric.J[ix,iy,iz,:,:,:],Cache)
      end
    end
  end  
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) 
  DSSC!(FC,JRhoC)
end

function AdvectionBlock!(Fc,c,uC,vC,Fe,dXdxI,J,Cache)
  DXW = Fe.DXW
  DYW = Fe.DYW
  @views W = Cache.Block[:,:,1,1]
  @views temp = Cache.Block[:,:,1,2]
  @views DxW = Cache.Block[:,:,1,2]
  @views DyW = Cache.Block[:,:,1,3]

  @views @. temp = 0.5 * c * ((dXdxI[:,:,1,1,1] + dXdxI[:,:,2,1,1]) * uC +
    (dXdxI[:,:,1,1,2] + dXdxI[:,:,2,1,2]) * vC)
  DerivativeX!(Fc,temp,DXW)
  @views @. temp = 0.5 * c * ((dXdxI[:,:,1,2,1] + dXdxI[:,:,2,2,1]) * uC +
    (dXdxI[:,:,1,2,2] + dXdxI[:,:,2,2,2]) * vC)
  DerivativeY!(Fc,temp,DYW)
end  

