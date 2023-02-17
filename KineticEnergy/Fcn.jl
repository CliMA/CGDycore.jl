function Fcn!(FC,FF,UC,UF,Metric,Fe,PhysParam,Cache)
  RhoCPos = 1
  uCPos = 2
  ThCPos = 3
  wFPos = 1

  Grav = PhysParam.Grav

  Nx = size(FC,1)
  Nz = size(FC,2)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
  pC = Cache.pC
  pF = Cache.pF
  uF = Cache.uF
  RhoF = Cache.RhoF
  RhoThetaF = Cache.RhoThetaF
  FuF = Cache.FuF
  FRhoF = Cache.FRhoF
  FRhoThetaF = Cache.FRhoThetaF
  KinF = Cache.KinF
  KinC = Cache.KinC

  @views wF = UF[:,:,:,:,wFPos]
  @views FwF = FF[:,:,:,:,wFPos]

  @views Pressure!(pC,UC[:,:,:,:,ThCPos],PhysParam)
  for ix = 1 : Nx
    for iz = 1 : Nz
      for i = 1 : OrdPolyX +1
        @views uF[ix,iz,i,:] = Fe.IntZC2F * UC[ix,iz,i,:,uCPos]
        @views RhoF[ix,iz,i,:] = Fe.IntZC2F * UC[ix,iz,i,:,RhoCPos]
        @views RhoThetaF[ix,iz,i,:] = Fe.IntZC2F * UC[ix,iz,i,:,ThCPos]
        @views pF[ix,iz,i,:] = Fe.IntZC2F * pC[ix,iz,i,:]
      end
    end
  end
  @. KinF = 0.5 * (wF * wF + uF * uF)
  for ix = 1 : Nx
    for iz = 1 : Nz
      for i = 1 : OrdPolyX +1
        @views KinC[ix,iz,i,:] = Fe.IntZF2C * (KinF[ix,iz,i,:] .+ Grav*Metric.X[ix,iz,i,:,2])
        @views KinF[ix,iz,i,:] = Fe.IntZC2F * KinC[ix,iz,i,:]
      end
    end
  end

  @. FuF = 0.0
  @. FwF = 0.0
  @. FRhoF = 0.0
  @. FRhoThetaF = 0.0
  Curl!(FuF,FwF,uF,wF,RhoF,Fe,Metric,Cache)
  Div!(FRhoF,uF,wF,RhoF,Fe,Metric)
  Div!(FRhoThetaF,uF,wF,RhoThetaF,Fe,Metric)
  RhoGrad!(FuF,FwF,KinF,RhoF,Fe,Metric)
  Grad!(FuF,FwF,pF,Fe,Metric)

  DSSF!(FwF,RhoF,Metric.J)
  DSS!(FuF,RhoF,Metric.J)
  DSS!(FRhoF,Metric.J)
  DSS!(FRhoThetaF,Metric.J)
  for ix = 1 : Nx
    for iz = 1 : Nz
      for i = 1 : OrdPolyX +1
        @views FC[ix,iz,i,:,uCPos] =  Fe.P * FuF[ix,iz,i,:]
        @views FC[ix,iz,i,:,RhoCPos] = Fe.P * FRhoF[ix,iz,i,:]
        @views FC[ix,iz,i,:,ThCPos] = Fe.P * FRhoThetaF[ix,iz,i,:]
      end
    end
  end
# Buoyancy!(FwF,PhysParam)
end  
