function DivRhoColumn!(FRhoC,uC,vC,w,RhoC,Fe,dXdxI,Cache)
  Nz = size(uC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views tempZ = Cache.Column[:,:,:,:,1]
  @views temp = Cache.Block[:,:,1,1]
  @views temp1 = Cache.Block[:,:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,1]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz  
    @views @. tempZ[:,:,1,iz] = -RhoC[iz,:,:] * (dXdxI[iz,:,:,1,3,1] * uC[iz,:,:] +
      +dXdxI[iz,:,:,1,3,2] * vC[iz,:,:] + dXdxI[iz,:,:,1,3,3] * w[iz,:,:])
    @views @. tempZ[:,:,2,iz] = -RhoC[iz,:,:] * (dXdxI[iz,:,:,2,3,1] * uC[iz,:,:] +
      +dXdxI[iz,:,:,2,3,2] * vC[iz,:,:] + dXdxI[iz,:,:,2,3,3] * w[iz+1,:,:])
  end    
  @inbounds for iz = 1 :Nz  
    @views @. temp = -RhoC[iz,:,:] * ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * uC[iz,:,:] +
      (dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,1,2]) * vC[iz,:,:] + 
      dXdxI[iz,:,:,1,1,3] * w[iz,:,:] + dXdxI[iz,:,:,2,1,3] * w[iz+1,:,:])
    @views DerivativeX!(FRhoC[iz,:,:],temp,DX)
    @views @. temp = -RhoC[iz,:,:] * ((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * uC[iz,:,:] +
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * vC[iz,:,:] + 
      dXdxI[iz,:,:,1,2,3] * w[iz,:,:] + dXdxI[iz,:,:,2,2,3] * w[iz+1,:,:])
    @views DerivativeY!(FRhoC[iz,:,:],temp,DY)
    @views @. temp1 = 0.0 
    @views DerivativeZ!(temp1,tempZ[:,:,:,iz],DZ)
    @views @. FRhoC[iz,:,:] += (temp1[:,:,1] + temp1[:,:,2]) 
    # Density
    if iz > 1
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz] - tempZ[:,:,2,iz-1]) 
      @views @. FRhoC[iz,:,:] += FluxZ
    end  
    if iz < Nz
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz+1] - tempZ[:,:,2,iz]) 
      @views @. FRhoC[iz,:,:] += FluxZ
    end  
  end 
end 

function DivRhoThetaColumn!(FRhoThetaC,uC,vC,w,RhoThetaC,Fe,dXdxI,Cache)
  Nz = size(uC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views tempZ = Cache.Column[:,:,:,:,1]
  @views temp = Cache.Block[:,:,1,1]
  @views temp1 = Cache.Block[:,:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,1]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz  
    @views @. tempZ[:,:,1,iz] = -RhoThetaC[iz,:,:] * (dXdxI[iz,:,:,1,3,1] * uC[iz,:,:] +
      +dXdxI[iz,:,:,1,3,2] * vC[iz,:,:] + dXdxI[iz,:,:,1,3,3] * w[iz,:,:])
    @views @. tempZ[:,:,2,iz] = -RhoThetaC[iz,:,:] * (dXdxI[iz,:,:,2,3,1] * uC[iz,:,:] +
      +dXdxI[iz,:,:,2,3,2] * vC[iz,:,:] + dXdxI[iz,:,:,2,3,3] * w[iz+1,:,:])
  end    
  @inbounds for iz = 1 :Nz  
    @views @. temp = -RhoThetaC[iz,:,:] * ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * uC[iz,:,:] +
      (dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,1,2]) * vC[iz,:,:] + 
      dXdxI[iz,:,:,1,1,3] * w[iz,:,:] + dXdxI[iz,:,:,2,1,3] * w[iz+1,:,:])
    @views DerivativeX!(FRhoThetaC[iz,:,:],temp,DX)
    @views @. temp = -RhoThetaC[iz,:,:] * ((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * uC[iz,:,:] +
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * vC[iz,:,:] + 
      dXdxI[iz,:,:,1,2,3] * w[iz,:,:] + dXdxI[iz,:,:,2,2,3] * w[iz+1,:,:])
    @views DerivativeY!(FRhoThetaC[iz,:,:],temp,DY)
    @views @. temp1 = 0.0 
    @views DerivativeZ!(temp1,tempZ[:,:,:,iz],DZ)
    @views @. FRhoThetaC[iz,:,:] += (temp1[:,:,1] + temp1[:,:,2]) 
    # Density
    if iz > 1
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz] - tempZ[:,:,2,iz-1]) 
      @views @. FRhoThetaC[iz,:,:] += FluxZ
    end  
    if iz < Nz
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz+1] - tempZ[:,:,2,iz]) 
      @views @. FRhoThetaC[iz,:,:] += FluxZ
    end  
  end 
end 

