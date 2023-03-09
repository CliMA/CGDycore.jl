function CurlColumn!(FuF,FvF,Fw,uF,vF,w,RhoF,Fe,dXdxI,Cache)
  Nz = size(w,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views tempuZ = Cache.Column[:,:,:,:,1]
  @views tempvZ = Cache.Column[:,:,:,:,2]
  @views tempwZ = Cache.Column[:,:,:,:,3]
  @views U = Cache.Block[:,:,:,1]
  @views V = Cache.Block[:,:,:,2]
  @views W = Cache.Block[:,:,:,3]
  @views temp = Cache.Block[:,:,:,4]

  @views FluxUZ = Cache.BlockXY[:,:,1]
  @views FluxVZ = Cache.BlockXY[:,:,2]
  @views FluxWZ = Cache.BlockXY[:,:,3]
  OPz = OrdPolyZ + 1
  @inbounds for iz = 1 : Nz 
     @views @. tempuZ[:,:,:,iz] = -dXdxI[iz,:,:,:,3,3] * vF[iz,:,:,:] + dXdxI[iz,:,:,:,3,2] * w[iz,:,:,:]
     @views @. tempvZ[:,:,:,iz] = dXdxI[iz,:,:,:,3,3] * uF[iz,:,:,:] - dXdxI[iz,:,:,:,3,1] * w[iz,:,:,:]
     @views @. tempwZ[:,:,:,iz] = dXdxI[iz,:,:,:,3,1] * vF[iz,:,:,:] - dXdxI[iz,:,:,:,3,2] * uF[iz,:,:,:]
  end  
  @inbounds for iz = 1 : Nz 
    if iz > 1
      @views @. FluxUZ = 0.5 * (tempuZ[:,:,1,iz] - tempuZ[:,:,OPz,iz-1])
      @views @. FluxVZ = 0.5 * (tempvZ[:,:,1,iz] - tempvZ[:,:,OPz,iz-1])
      @views @. FluxWZ = 0.5 * (tempwZ[:,:,1,iz] - tempwZ[:,:,OPz,iz-1])
      @views @. FuF[iz,:,:,1] += RhoF[iz,:,:,1] * (vF[iz,:,:,1] * FluxWZ - w[iz,:,:,1] * FluxVZ)
      @views @. FvF[iz,:,:,1] += RhoF[iz,:,:,1] * (uF[iz,:,:,1] * FluxWZ - w[iz,:,:,1] * FluxUZ)
      @views @. Fw[iz,:,:,1] += RhoF[iz,:,:,1] * (uF[iz,:,:,1] * FluxVZ - vF[iz,:,:,1] * FluxUZ)
    end
    if iz < Nz
      @views @. FluxUZ = 0.5 * (tempuZ[:,:,1,iz+1] - tempuZ[:,:,OPz,iz])
      @views @. FluxVZ = 0.5 * (tempvZ[:,:,1,iz+1] - tempvZ[:,:,OPz,iz])
      @views @. FluxWZ = 0.5 * (tempwZ[:,:,1,iz+1] - tempwZ[:,:,OPz,iz])
      @views @. FuF[iz,:,:,OPz] += RhoF[iz,:,:,OPz] * (vF[iz,:,:,OPz] * FluxWZ - w[iz,:,:,OPz] * FluxVZ)
      @views @. FvF[iz,:,:,OPz] += RhoF[iz,:,:,OPz] * (uF[iz,:,:,OPz] * FluxWZ - w[iz,:,:,OPz] * FluxUZ)
      @views @. Fw[iz,:,:,OPz] += RhoF[iz,:,:,OPz] * (uF[iz,:,:,OPz] * FluxVZ - vF[iz,:,:,OPz] * FluxUZ)
    end    
  
    @. U = 0.0
    @. V = 0.0
    @. W = 0.0
    @views @. temp = -dXdxI[iz,:,:,:,1,3] * vF[iz,:,:,:] + dXdxI[iz,:,:,:,1,2] * w[iz,:,:,:]
    DerivativeX!(U,temp,DX)
    @views @. temp = dXdxI[iz,:,:,:,1,3] * uF[iz,:,:,:] - dXdxI[iz,:,:,:,1,1] * w[iz,:,:,:]
    DerivativeX!(V,temp,DX)
    @views @. temp = dXdxI[iz,:,:,:,1,1] * vF[iz,:,:,:] - dXdxI[iz,:,:,:,1,2] * uF[iz,:,:,:]
    DerivativeX!(W,temp,DX)
    @views @. temp = -dXdxI[iz,:,:,:,2,3] * vF[iz,:,:,:] + dXdxI[iz,:,:,:,2,2] * w[iz,:,:,:]
    DerivativeY!(U,temp,DY)
    @views @. temp = dXdxI[iz,:,:,:,2,3] * uF[iz,:,:,:] - dXdxI[iz,:,:,:,2,1] * w[iz,:,:,:]
    DerivativeY!(V,temp,DY)
    @views @. temp = dXdxI[iz,:,:,:,2,1] * vF[iz,:,:,:] - dXdxI[iz,:,:,:,2,2] * uF[iz,:,:,:]
    DerivativeY!(W,temp,DY)
    @views DerivativeZ!(U,tempuZ[:,:,:,iz],DZ)
    @views DerivativeZ!(V,tempvZ[:,:,:,iz],DZ)
    @views DerivativeZ!(W,tempwZ[:,:,:,iz],DZ)
    @views @. FuF[iz,:,:,:] += RhoF[iz,:,:,:] * (vF[iz,:,:,:] * W - w[iz,:,:,:] * V)
    @views @. FvF[iz,:,:,:] += RhoF[iz,:,:,:] * (uF[iz,:,:,:] * W - w[iz,:,:,:] * U)
    @views @. Fw[iz,:,:,:] += RhoF[iz,:,:,:] * (uF[iz,:,:,:] * V - vF[iz,:,:,:] * U)
      
  end 
end

