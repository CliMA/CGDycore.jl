function DivRhoColumn!(FRhoF,uF,vF,w,RhoF,Fe,dXdxI,Cache)
  Nz = size(w,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views tempZ = Cache.Column[:,:,:,:,1]
  @views temp = Cache.Block[:,:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,1]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz  
    @views @. tempZ[:,:,:,iz] = -RhoF[iz,:,:,:] * (dXdxI[iz,:,:,:,3,1] * uF[iz,:,:,:] +
      +dXdxI[iz,:,:,:,3,2] * vF[iz,:,:,:] + dXdxI[iz,:,:,:,3,3] * w[iz,:,:,:])
  end    
  @inbounds for iz = 1 :Nz  
    @views @. temp = -RhoF[iz,:,:,:] * (dXdxI[iz,:,:,:,1,1] * uF[iz,:,:,:] +
      dXdxI[iz,:,:,:,1,2] * vF[iz,:,:,:] + dXdxI[iz,:,:,:,1,3] * w[iz,:,:,:])
    @views DerivativeX!(FRhoF[iz,:,:,:],temp,DX)
    @views @. temp = -RhoF[iz,:,:,:] * (dXdxI[iz,:,:,:,2,1] * uF[iz,:,:,:] +
      dXdxI[iz,:,:,:,2,2] * vF[iz,:,:,:] + dXdxI[iz,:,:,:,2,3] * w[iz,:,:,:])
    @views DerivativeY!(FRhoF[iz,:,:,:],temp,DY)
    @views DerivativeZ!(FRhoF[iz,:,:,:],tempZ[:,:,:,iz],DZ)
    # Density
    if iz > 1
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz] - tempZ[:,:,OPz,iz-1]) 
      @views @. FRhoF[iz,:,:,1] += FluxZ
    end  
    if iz < Nz
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz+1] -tempZ[:,:,OPz,iz]) 
      @views @. FRhoF[iz,:,:,OPz] += FluxZ
    end  
  end 
end 

function DivRhoThetaColumn!(FRhoThetaF,uF,vF,w,RhoThetaF,Fe,dXdxI,Cache)
  Nz = size(w,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views tempZ = Cache.Column[:,:,:,:,1]
  @views temp = Cache.Block[:,:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,1]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz  
    @views @. tempZ[:,:,:,iz] = -RhoThetaF[iz,:,:,:] * (dXdxI[iz,:,:,:,3,1] * uF[iz,:,:,:] +
      +dXdxI[iz,:,:,:,3,2] * vF[iz,:,:,:] + dXdxI[iz,:,:,:,3,3] * w[iz,:,:,:])
  end    
  @inbounds for iz = 1 : Nz  
    @views @. temp = -RhoThetaF[iz,:,:,:] * (dXdxI[iz,:,:,:,1,1] * uF[iz,:,:,:] +
      dXdxI[iz,:,:,:,1,2] * vF[iz,:,:,:] + dXdxI[iz,:,:,:,1,3] * w[iz,:,:,:])
    @views DerivativeX!(FRhoThetaF[iz,:,:,:],temp,DX)
    @views @. temp = -RhoThetaF[iz,:,:,:] * (dXdxI[iz,:,:,:,2,1] * uF[iz,:,:,:] +
      dXdxI[iz,:,:,:,2,2] * vF[iz,:,:,:] + dXdxI[iz,:,:,:,2,3] * w[iz,:,:,:])
    @views DerivativeY!(FRhoThetaF[iz,:,:,:],temp,DY)
    @views DerivativeZ!(FRhoThetaF[iz,:,:,:],tempZ[:,:,:,iz],DZ)

    # Density
    if iz > 1
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz] - tempZ[:,:,OPz,iz-1]) 
      @views @. FRhoThetaF[iz,:,:,1] += FluxZ
    end  
    if iz < Nz
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz+1] -tempZ[:,:,OPz,iz]) 
      @views @. FRhoThetaF[iz,:,:,OPz] += FluxZ
    end  
  end 
end 

