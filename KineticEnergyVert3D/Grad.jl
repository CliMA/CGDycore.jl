function RhoGradColumn!(FuF,FvF,Fw,pF,RhoF,Fe,dXdxI,Cache)
  Nz = size(pF,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views DXpF = Cache.Block[:,:,:,1]
  @views DYpF = Cache.Block[:,:,:,2]
  @views DZpF = Cache.Block[:,:,:,3]
  @views GradZ = Cache.BlockXY[:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,2]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz 
    @. DXpF = 0.0
    @views DerivativeX!(DXpF,pF[iz,:,:,:],DX)
    @. DYpF = 0.0
    @views DerivativeY!(DYpF,pF[iz,:,:,:],DY)
    @. DZpF = 0.0
    @views DerivativeZ!(DZpF,pF[iz,:,:,:],DZ)
    @views @. FuF[iz,:,:,:] -= RhoF[iz,:,:,:] * 
      (dXdxI[iz,:,:,:,1,1] * DXpF + dXdxI[iz,:,:,:,2,1] * DYpF + dXdxI[iz,:,:,:,3,1] * DZpF[:,:,:])
    @views @. FvF[iz,:,:,:] -= RhoF[iz,:,:,:] * 
      (dXdxI[iz,:,:,:,1,2] * DXpF + dXdxI[iz,:,:,:,2,2] * DYpF + dXdxI[iz,:,:,:,3,2] * DZpF[:,:,:])
    @views @. Fw[iz,:,:,:] -= RhoF[iz,:,:,:] * 
      (dXdxI[iz,:,:,:,1,3] * DXpF + dXdxI[iz,:,:,:,2,3] * DYpF + dXdxI[iz,:,:,:,3,3] * DZpF[:,:,:])
    if iz > 1
      @views @. GradZ = 0.5 * (pF[iz,:,:,1] - pF[iz-1,:,:,OPz]) * RhoF[iz,:,:,1]  
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,1]
      @views @. FuF[iz,:,:,1] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,2]
      @views @. FvF[iz,:,:,1] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,3]
      @views @. Fw[iz,:,:,1] -= FluxZ 
    end    
    if iz < Nz 
      @views @. GradZ = 0.5 * (pF[iz+1,:,:,1] - pF[iz,:,:,OPz]) * RhoF[iz,:,:,OPz] 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,OPz,3,1] 
      @views @. FuF[iz,:,:,OPz] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,OPz,3,2]
      @views @. FvF[iz,:,:,OPz] -= FluxZ 
      @views @. FluxZ =  GradZ * dXdxI[iz,:,:,OPz,3,3]
      @views @. Fw[iz,:,:,OPz] -= FluxZ 
    end    
  end 
end 
function GradColumn!(FuF,FvF,Fw,pF,Fe,dXdxI,Cache)
  Nz = size(pF,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views DXpF = Cache.Block[:,:,:,1]
  @views DYpF = Cache.Block[:,:,:,2]
  @views DZpF = Cache.Block[:,:,:,3]
  @views GradZ = Cache.BlockXY[:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,2]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz 
    @. DXpF = 0.0
    @views DerivativeX!(DXpF,pF[iz,:,:,:],DX)
    @. DYpF = 0.0
    @views DerivativeY!(DYpF,pF[iz,:,:,:],DY)
    @. DZpF = 0.0
    @views DerivativeZ!(DZpF,pF[iz,:,:,:],DZ)
    @views @. FuF[iz,:,:,:] -= 
      (dXdxI[iz,:,:,:,1,1] * DXpF + dXdxI[iz,:,:,:,2,1] * DYpF + dXdxI[iz,:,:,:,3,1] * DZpF[:,:,:])
    @views @. FvF[iz,:,:,:] -= 
      (dXdxI[iz,:,:,:,1,2] * DXpF + dXdxI[iz,:,:,:,2,2] * DYpF + dXdxI[iz,:,:,:,3,2] * DZpF[:,:,:])
    @views @. Fw[iz,:,:,:] -= 
      (dXdxI[iz,:,:,:,1,3] * DXpF + dXdxI[iz,:,:,:,2,3] * DYpF + dXdxI[iz,:,:,:,3,3] * DZpF[:,:,:])
    if iz > 1
      @views @. GradZ = 0.5 * (pF[iz,:,:,1] - pF[iz-1,:,:,OPz])  
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,1]
      @views @. FuF[iz,:,:,1] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,2]
      @views @. FvF[iz,:,:,1] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,1,3,3]
      @views @. Fw[iz,:,:,1] -= FluxZ 
    end    
    if iz < Nz 
      @views @. GradZ = 0.5 * (pF[iz+1,:,:,1] - pF[iz,:,:,OPz])  
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,OPz,3,1]
      @views @. FuF[iz,:,:,OPz] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[iz,:,:,OPz,3,2]
      @views @. FvF[iz,:,:,OPz] -= FluxZ 
      @views @. FluxZ =  GradZ * dXdxI[iz,:,:,OPz,3,3]
      @views @. Fw[iz,:,:,OPz] -= FluxZ 
    end    
  end 
end 

