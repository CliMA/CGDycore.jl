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
  @views uF = Cache.Block[:,:,1,2]
  @views vF = Cache.Block[:,:,1,3]
  @views JXY = Cache.CF[:,:,1,:,:,1,1]
  for ix = 1 : Nx
    for iy = 1 : Ny
      for i = 1 : OrdPolyX + 1  
        for j = 1 : OrdPolyY + 1  
          uF[i,j] = 0.0  
          vF[i,j] = 0.0  
          for l = 1 : OrdPolyZ   
            uF[i,j] += Fe.IntZC2F[1,l] * UC[ix,iy,1,i,j,l,uPos]
            vF[i,j] += Fe.IntZC2F[1,l] * UC[ix,iy,1,i,j,l,vPos]
          end  
        end
      end  
      @views @. wF[ix,iy,1,:,:,1] = -(dXdxI[ix,iy,1,:,:,1,3,1]*uF[:,:] +
        dXdxI[ix,iy,1,:,:,1,3,2]*vF[:,:]) /
        dXdxI[ix,iy,1,:,:,1,3,3]  
    end    
  end    
  #Average wF
  @. JXY = 1.0
  @views AverageFXY!(wF[:,:,1,:,:,1],JXY)
end
