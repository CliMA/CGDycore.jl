function RhoGradColumn!(FuC,FvC,Fw,pC,RhoC,Fe,dXdxI,Cache)
  Nz = size(pC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views DXpC = Cache.Block[:,:,1,1]
  @views DYpC = Cache.Block[:,:,1,2]
  @views DZpC = Cache.Block[:,:,1,3]
  @views GradZ = Cache.BlockXY[:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,2]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz 
    @. DXpC = 0.0
    @views DerivativeX!(DXpC,pC[iz,:,:],DX)
    @. DYpC = 0.0
    @views DerivativeY!(DYpC,pC[iz,:,:],DY)
#   @. DZpC = 0.0
#   @views DerivativeZ!(DZpC,pC[iz,:,:],DZ)
    @views @. FuC[iz,:,:] -= 0.5 * RhoC[iz,:,:] * 
      ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * DXpC + 
      (dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,2,1]) * DYpC)
    @views @. FvC[iz,:,:] -=  0.5 * RhoC[iz,:,:] *
      ((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * DXpC + 
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * DYpC)
    @views @. Fw[iz,:,:] -= RhoC[iz,:,:] * (dXdxI[iz,:,:,1,1,3] * DXpC + dXdxI[iz,:,:,1,2,3] * DYpC)
    @views @. Fw[iz+1,:,:] -= RhoC[iz,:,:] * (dXdxI[iz,:,:,2,1,3] * DXpC + dXdxI[iz,:,:,2,2,3] * DYpC)
    if iz > 1
      @views @. GradZ = 0.5 * (pC[iz,:,:] - pC[iz-1,:,:]) * RhoC[iz,:,:] 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,1]
      @views @. FuC[iz,:,:] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,2]
      @views @. FvC[iz,:,:] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,3]
      @views @. Fw[iz,:,:] -= 0.5 * FluxZ 
    end    
    if iz < Nz 
      @views @. GradZ = 0.5 * (pC[iz+1,:,:] - pC[iz,:,:]) * RhoC[iz,:,:]  
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,2,3,1]
      @views @. FuC[iz,:,:] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,2,3,2]
      @views @. FvC[iz,:,:] -= FluxZ 
      @views @. FluxZ =  GradZ * dXdxI[iz,:,:,2,3,3]
      @views @. Fw[iz+1,:,:] -= 0.5 * FluxZ 
    end    
  end 
end 

function GradColumn!(FuC,FvC,Fw,pC,Fe,dXdxI,Cache)
  Nz = size(pC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views DXpC = Cache.Block[:,:,1,1]
  @views DYpC = Cache.Block[:,:,1,2]
  @views DZpC = Cache.Block[:,:,1,3]
  @views GradZ = Cache.BlockXY[:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,2]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz 
    @. DXpC = 0.0
    @views DerivativeX!(DXpC,pC[iz,:,:],DX)
    @. DYpC = 0.0
    @views DerivativeY!(DYpC,pC[iz,:,:],DY)
    @views @. FuC[iz,:,:] -= 0.5 * 
      ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * DXpC + 
      (dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,2,1]) * DYpC)
    @views @. FvC[iz,:,:] -=  0.5 * 
      ((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * DXpC + 
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * DYpC)
    @views @. Fw[iz,:,:] -= (dXdxI[iz,:,:,1,1,3] * DXpC + dXdxI[iz,:,:,1,2,3] * DYpC)
    @views @. Fw[iz+1,:,:] -= (dXdxI[iz,:,:,2,1,3] * DXpC + dXdxI[iz,:,:,2,2,3] * DYpC)
    if iz > 1
      @views @. GradZ = 0.5 * (pC[iz,:,:] - pC[iz-1,:,:])
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,1]
      @views @. FuC[iz,:,:] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,2]
      @views @. FvC[iz,:,:] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,3]
      @views @. Fw[iz,:,:] -= 0.5 * FluxZ 
    end    
    if iz < Nz 
      @views @. GradZ = 0.5 * (pC[iz+1,:,:] - pC[iz,:,:])
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,2,3,1]
      @views @. FuC[iz,:,:] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,2,3,2]
      @views @. FvC[iz,:,:] -= FluxZ 
      @views @. FluxZ =  GradZ * dXdxI[iz,:,:,2,3,3]
      @views @. Fw[iz+1,:,:] -= 0.5 * FluxZ 
    end    
  end 
end 

