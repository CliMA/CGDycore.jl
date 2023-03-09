function CurlColumn!(FuC,FvC,Fw,uC,vC,w,RhoC,Fe,dXdxI,Cache)
  Nz = size(FuC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views tempuZ = Cache.Column[:,:,1:2,:,1]
  @views tempvZ = Cache.Column[:,:,1:2,:,2]
  @views tempwZ = Cache.Column[:,:,1:2,:,3]
  @views U = Cache.Block[:,:,:,1]
  @views V = Cache.Block[:,:,:,2]
  @views W = Cache.Block[:,:,:,3]
  @views temp = Cache.Block[:,:,1:2,4]

  @views FluxUZ = Cache.BlockXY[:,:,1]
  @views FluxVZ = Cache.BlockXY[:,:,2]
  @views FluxWZ = Cache.BlockXY[:,:,3]
  OPz = OrdPolyZ + 1
  @inbounds for iz = 1 : Nz 
    @views @. tempuZ[:,:,1,iz] = -dXdxI[iz,:,:,1,3,3] * vC[iz,:,:] + dXdxI[iz,:,:,1,3,2] * w[iz,:,:]
    @views @. tempuZ[:,:,2,iz] = -dXdxI[iz,:,:,2,3,3] * vC[iz,:,:] + dXdxI[iz,:,:,2,3,2] * w[iz+1,:,:]
    @views @. tempvZ[:,:,1,iz] = dXdxI[iz,:,:,1,3,3] * uC[iz,:,:] - dXdxI[iz,:,:,1,3,1] * w[iz,:,:]
    @views @. tempvZ[:,:,2,iz] = dXdxI[iz,:,:,2,3,3] * uC[iz,:,:] - dXdxI[iz,:,:,2,3,1] * w[iz+1,:,:]
    @views @. tempwZ[:,:,1,iz] = dXdxI[iz,:,:,1,3,1] * vC[iz,:,:] - dXdxI[iz,:,:,1,3,2] * uC[iz,:,:]
    @views @. tempwZ[:,:,2,iz] = dXdxI[iz,:,:,2,3,1] * vC[iz,:,:] - dXdxI[iz,:,:,2,3,2] * uC[iz,:,:]
  end  
  @inbounds for iz = 1 : Nz 
    if iz > 1
      @views @. FluxUZ = 0.5 * (tempuZ[:,:,1,iz] - tempuZ[:,:,2,iz-1])
      @views @. FluxVZ = 0.5 * (tempvZ[:,:,1,iz] - tempvZ[:,:,2,iz-1])
      @views @. FluxWZ = 0.5 * (tempwZ[:,:,1,iz] - tempwZ[:,:,2,iz-1])
      @views @. FuC[iz,:,:] += RhoC[iz,:,:] * (vC[iz,:,:,1] * FluxWZ - w[iz,:,:,1] * FluxVZ)
      @views @. FvC[iz,:,:] += RhoC[iz,:,:] * (uC[iz,:,:,1] * FluxWZ - w[iz,:,:,1] * FluxUZ)
      @views @. Fw[iz,:,:] += RhoC[iz,:,:] * (uC[iz,:,:,1] * FluxVZ - vC[iz,:,:,1] * FluxUZ)
    end
    if iz < Nz
      @views @. FluxUZ = 0.5 * (tempuZ[:,:,1,iz+1] - tempuZ[:,:,2,iz])
      @views @. FluxVZ = 0.5 * (tempvZ[:,:,1,iz+1] - tempvZ[:,:,2,iz])
      @views @. FluxWZ = 0.5 * (tempwZ[:,:,1,iz+1] - tempwZ[:,:,2,iz])
      @views @. FuC[iz,:,:] += RhoC[iz,:,:] * (vC[iz,:,:,] * FluxWZ - w[iz+1,:,:] * FluxVZ)
      @views @. FvC[iz,:,:] += RhoC[iz,:,:] * (uC[iz,:,:,] * FluxWZ - w[iz+1,:,:] * FluxUZ)
      @views @. Fw[iz+1,:,:] += RhoC[iz,:,:,] * (uC[iz,:,:] * FluxVZ - vC[iz,:,:] * FluxUZ)
    end    
  
    @. U = 0.0
    @. V = 0.0
    @. W = 0.0
    @views @. temp[:,:,1] = -dXdxI[iz,:,:,1,1,3] * vC[iz,:,:] + dXdxI[iz,:,:,1,1,2] * w[iz,:,:]
    @views @. temp[:,:,2] = -dXdxI[iz,:,:,2,1,3] * vC[iz,:,:] + dXdxI[iz,:,:,2,1,2] * w[iz+1,:,:]
    DerivativeX!(U,temp,DX)
    @views @. temp[:,:,1] = dXdxI[iz,:,:,1,1,3] * uC[iz,:,:] - dXdxI[iz,:,:,1,1,1] * w[iz,:,:]
    @views @. temp[:,:,2] = dXdxI[iz,:,:,2,1,3] * uC[iz,:,:] - dXdxI[iz,:,:,2,1,1] * w[iz+1,:,:]
    DerivativeX!(V,temp,DX)
    @views @. temp[:,:,1] = dXdxI[iz,:,:,1,1,1] * vC[iz,:,:] - dXdxI[iz,:,:,1,1,2] * uC[iz,:,:]
    @views @. temp[:,:,2] = dXdxI[iz,:,:,2,1,1] * vC[iz,:,:] - dXdxI[iz,:,:,2,1,2] * uC[iz,:,:]
    DerivativeX!(W,temp,DX)
    @views @. temp[:,:,1] = -dXdxI[iz,:,:,1,2,3] * vC[iz,:,:] + dXdxI[iz,:,:,1,2,2] * w[iz,:,:]
    @views @. temp[:,:,2] = -dXdxI[iz,:,:,2,2,3] * vC[iz,:,:] + dXdxI[iz,:,:,2,2,2] * w[iz+1,:,:]
    DerivativeY!(U,temp,DY)
    @views @. temp[:,:,1] = dXdxI[iz,:,:,1,2,3] * uC[iz,:,:] - dXdxI[iz,:,:,1,2,1] * w[iz,:,:]
    @views @. temp[:,:,2] = dXdxI[iz,:,:,2,2,3] * uC[iz,:,:] - dXdxI[iz,:,:,2,2,1] * w[iz+1,:,:]
    DerivativeY!(V,temp,DY)
    @views @. temp[:,:,1] = dXdxI[iz,:,:,1,2,1] * vC[iz,:,:] - dXdxI[iz,:,:,1,2,2] * uC[iz,:,:]
    @views @. temp[:,:,2] = dXdxI[iz,:,:,2,2,1] * vC[iz,:,:] - dXdxI[iz,:,:,2,2,2] * uC[iz,:,:]
    DerivativeY!(W,temp,DY)
    @views DerivativeZ!(U,tempuZ[:,:,:,iz],DZ)
    @views DerivativeZ!(V,tempvZ[:,:,:,iz],DZ)
    @views DerivativeZ!(W,tempwZ[:,:,:,iz],DZ)
    @views @. FuC[iz,:,:] += RhoC[iz,:,:] * (vC[iz,:,:] * W[:,:,1] - w[iz,:,:] * V[:,:,1]) +
      RhoC[iz,:,:] * (vC[iz,:,:] * W[:,:,2] - w[iz+1,:,:] * V[:,:,2]) 
    @views @. FvC[iz,:,:] += RhoC[iz,:,:] * (uC[iz,:,:] * W[:,:,1] - w[iz,:,:,:] * U[:,:,1]) +
      RhoC[iz,:,:] * (uC[iz,:,:] * W[:,:,2] - w[iz+1,:,:,:] * U[:,:,2]) +
    @views @. Fw[iz,:,:] += RhoC[iz,:,:] * (uC[iz,:,:] * V[:,:,1] - vC[iz,:,:] * U[:,:,1])
    @views @. Fw[iz+1,:,:] += RhoC[iz,:,:] * (uC[iz,:,:] * V[:,:,2] - vC[iz,:,:] * U[:,:,2])
      
  end 
end

