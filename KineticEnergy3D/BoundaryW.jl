function BoundaryW!(wF,UC,dXdxI,Fe,Cache)  
# with Interpolation for higher order
  uPos = 2
  vPos = 3
  Nx = size(wF,1)
  Ny = size(wF,2)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
# Lower boundary condition
  @views uC = UC[:,:,1,:,:,2]
  @views vC = UC[:,:,1,:,:,3]
  @views JXY = Cache.CF[:,:,1,:,:,1,1]
  for ix = 1 : Nx
    for iy = 1 : Ny
      @views @. wF[ix,iy,1,:,:] = -(dXdxI[ix,iy,1,:,:,1,3,1]*uC[ix,iy,:,:] +
        dXdxI[ix,iy,1,:,:,1,3,2]*vC[ix,iy,:,:]) /
        dXdxI[ix,iy,1,:,:,1,3,3]  
    end    
  end    
  #Average wF
  #@. JXY = 1.0
  #@views AverageFXY!(wF[:,:,1,:,:],JXY)
end
