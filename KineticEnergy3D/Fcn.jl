function Curl!(FC,FF,UC,UF,Metric,Fe,PhysParam,Cache)
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
  @views RhoF = Cache.UUF[:,:,:,:,:,:,1]
  @views uF = Cache.UUF[:,:,:,:,:,:,2]
  @views vF = Cache.UUF[:,:,:,:,:,:,3]
  @views FuF = Cache.FUUF[:,:,:,:,:,:,2]
  @views FvF = Cache.FUUF[:,:,:,:,:,:,3]

  @views wF = UF[:,:,:,:,:,:,wFPos]
  @views FwF = FF[:,:,:,:,:,:,wFPos]

  JRhoF = Cache.KinF

  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @inbounds for iz = 1 : Nz
        @inbounds for i = 1 : OrdPolyX +1
          @inbounds for j = 1 : OrdPolyY +1
            @views uF[ix,iy,iz,i,j,:] = Fe.IntZC2F * UC[ix,iy,iz,i,j,:,uCPos]
            @views RhoF[ix,iy,iz,i,j,:] = Fe.IntZC2F * UC[ix,iy,iz,i,j,:,RhoCPos]
          end
        end
      end
    end
  end
  @. FuF = 0.0
  @. FvF = 0.0
  @. FwF = 0.0
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @views CurlColumn!(FuF[ix,iy,:,:,:,:],FvF[ix,iy,:,:,:,:],FwF[ix,iy,:,:,:,:],
        uF[ix,iy,:,:,:,:],vF[ix,iy,:,:,:,:],wF[ix,iy,:,:,:,:],
        RhoF[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
    end
  end
  @. JRhoF = Metric.J * RhoF
  DSSF!(FwF,JRhoF)
  @. JRhoF = Metric.J * RhoF
  DSS!(FuF,JRhoF)
  @. JRhoF = Metric.J * RhoF
  DSS!(FvF,JRhoF)
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @inbounds for iz = 1 : Nz
        @inbounds for i = 1 : OrdPolyX +1
          @inbounds for j = 1 : OrdPolyY +1
            @views FC[ix,iy,iz,i,j,:,uCPos] =  Fe.P * FuF[ix,iy,iz,i,j,:]
            @views FC[ix,iy,iz,i,j,:,vCPos] =  Fe.P * FvF[ix,iy,iz,i,j,:]
          end
        end
      end
    end
  end
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
  @views RhoF = Cache.UUF[:,:,:,:,:,:,1]
  @views uF = Cache.UUF[:,:,:,:,:,:,2]
  @views vF = Cache.UUF[:,:,:,:,:,:,3]

  @views wF = UF[:,:,:,:,:,:,wFPos]

  @views FRhoF = Cache.FUUF[:,:,:,:,:,:,1]

  JRhoF = Cache.KinF

  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @inbounds for iz = 1 : Nz
        @inbounds for i = 1 : OrdPolyX +1
          @inbounds for j = 1 : OrdPolyY +1
            @views uF[ix,iy,iz,i,j,:] = Fe.IntZC2F * UC[ix,iy,iz,i,j,:,uCPos]
            @views vF[ix,iy,iz,i,j,:] = Fe.IntZC2F * UC[ix,iy,iz,i,j,:,vCPos]
            @views RhoF[ix,iy,iz,i,j,:] = Fe.IntZC2F * UC[ix,iy,iz,i,j,:,RhoCPos]
          end
        end
      end
    end
  end
  @. FRhoF = 0.0
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @views DivColumn!(FRhoF[ix,iy,:,:,:,:],uF[ix,iy,:,:,:,:],vF[ix,iy,:,:,:,:],
        wF[ix,iy,:,:,:,:],RhoF[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
    end
  end
  @. JRhoF = Metric.J 
  DSS!(FRhoF,JRhoF)
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @inbounds for iz = 1 : Nz
        @inbounds for i = 1 : OrdPolyX +1
          @inbounds for j = 1 : OrdPolyY +1
            @views FC[ix,iy,iz,i,j,:,RhoCPos] =  Fe.P * FRhoF[ix,iy,iz,i,j,:]
          end
        end
      end
    end
  end
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

  @views RhoF = Cache.UUF[:,:,:,:,:,:,1]
  @views uF = Cache.UUF[:,:,:,:,:,:,2]
  @views vF = Cache.UUF[:,:,:,:,:,:,3]
  @views wF = UF[:,:,:,:,:,:,wFPos]

  @views FuF = Cache.FUUF[:,:,:,:,:,:,2]
  @views FvF = Cache.FUUF[:,:,:,:,:,:,3]
  @views FwF = FF[:,:,:,:,:,:,wFPos]

  KinF = Cache.KinF
  JRhoF = Cache.KinF
  KinC = Cache.KinC

  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @inbounds for iz = 1 : Nz
        @inbounds for i = 1 : OrdPolyX +1
          @inbounds for j = 1 : OrdPolyY +1
            @views RhoF[ix,iy,iz,i,j,:] = Fe.IntZC2F * UC[ix,iy,iz,i,j,:,RhoCPos]
            @views uF[ix,iy,iz,i,j,:] = Fe.IntZC2F * UC[ix,iy,iz,i,j,:,uCPos]
            @views vF[ix,iy,iz,i,j,:] = Fe.IntZC2F * UC[ix,iy,iz,i,j,:,vCPos]
          end
        end
      end
    end
  end
  @. KinF = 0.5 * (wF * wF + vF * vF + uF * uF)
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @inbounds for iz = 1 : Nz
        @inbounds for i = 1 : OrdPolyX +1
          @inbounds for j = 1 : OrdPolyY +1
            @views KinC[ix,iy,iz,i,j,:] = Fe.IntZF2C * KinF[ix,iy,iz,i,j,:] 
            @views KinF[ix,iy,iz,i,j,:] = Fe.IntZC2F * KinC[ix,iy,iz,i,j,:]
          end
        end
      end
    end
  end
  @. FuF = 0.0
  @. FvF = 0.0
  @. FwF = 0.0
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @views RhoGradColumn!(FuF[ix,iy,:,:,:,:],FvF[ix,iy,:,:,:,:],FwF[ix,iy,:,:,:,:],
        KinF[ix,iy,:,:,:,:],RhoF[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
    end
  end
  @. JRhoF = Metric.J * RhoF
  DSSF!(FwF,JRhoF)
  @. JRhoF = Metric.J * RhoF
  DSS!(FuF,JRhoF)
  @. JRhoF = Metric.J * RhoF
  DSS!(FvF,JRhoF)
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @inbounds for iz = 1 : Nz
        @inbounds for i = 1 : OrdPolyX +1
          @inbounds for j = 1 : OrdPolyY +1
            @views FC[ix,iy,iz,i,j,:,uCPos] =  Fe.P * FuF[ix,iy,iz,i,j,:]
            @views FC[ix,iy,iz,i,j,:,vCPos] =  Fe.P * FvF[ix,iy,iz,i,j,:]
          end
        end
      end
    end
  end
end

function Fcn!(FC,FF,UC,UF,Metric,Fe,PhysParam,Cache)

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
  KinF = Cache.KinF
  JRhoF = Cache.KinF
  KinC = Cache.KinC

  @views wF = UF[:,:,:,:,:,:,wFPos]
  @views FwF = FF[:,:,:,:,:,:,wFPos]

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
  @. KinF = 0.5 * (wF * wF + vF * vF + uF * uF)
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @inbounds for iz = 1 : Nz
        @inbounds for i = 1 : OrdPolyX +1
          @inbounds for j = 1 : OrdPolyY +1
            @inbounds for k = 1 : OrdPolyZ
              KinC[ix,iy,iz,i,j,k] = 0.0
              @inbounds for l = 1 : OrdPolyZ + 1
                KinC[ix,iy,iz,i,j,k] += IntZF2C[k,l] * (KinF[ix,iy,iz,i,j,l] .+ Grav*Metric.X[ix,iy,iz,i,j,l,3])
              end  
            end  
            @inbounds for k = 1 : OrdPolyZ + 1 
              KinF[ix,iy,iz,i,j,k] = 0.0
              @inbounds for l = 1 : OrdPolyZ
                KinF[ix,iy,iz,i,j,k] += IntZC2F[k,l] * KinC[ix,iy,iz,i,j,l]
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
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @views CurlColumn!(FuF[ix,iy,:,:,:,:],FvF[ix,iy,:,:,:,:],FwF[ix,iy,:,:,:,:],
        uF[ix,iy,:,:,:,:],vF[ix,iy,:,:,:,:],wF[ix,iy,:,:,:,:],
        RhoF[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
      @views RhoGradColumn!(FuF[ix,iy,:,:,:,:],FvF[ix,iy,:,:,:,:],FwF[ix,iy,:,:,:,:],
        KinF[ix,iy,:,:,:,:],RhoF[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
      @views GradColumn!(FuF[ix,iy,:,:,:,:],FvF[ix,iy,:,:,:,:],FwF[ix,iy,:,:,:,:],
        pF[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
      @views DivColumn!(FRhoF[ix,iy,:,:,:,:],uF[ix,iy,:,:,:,:],vF[ix,iy,:,:,:,:],
        wF[ix,iy,:,:,:,:],RhoF[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
      @views DivColumn!(FRhoThetaF[ix,iy,:,:,:,:],uF[ix,iy,:,:,:,:],vF[ix,iy,:,:,:,:],
        wF[ix,iy,:,:,:,:],RhoThetaF[ix,iy,:,:,:,:],Fe,Metric.dXdxI[ix,iy,:,:,:,:,:,:],Cache)
    end
  end  

  @. JRhoF = Metric.J * RhoF
  DSSF!(FwF,JRhoF)
  @. JRhoF = Metric.J * RhoF
  DSS!(FuF,JRhoF)
  @. JRhoF = Metric.J * RhoF
  DSS!(FvF,JRhoF)
  @. JRhoF = Metric.J 
  DSS!(FRhoF,JRhoF)
  @. JRhoF = Metric.J 
  DSS!(FRhoThetaF,JRhoF)
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
# Buoyancy!(FwF,PhysParam)
end  
