function Curl!(FuLGL,Fw,uLGL,w,RhoLGL,Fe,Metric,Cache)
  Nx = size(w,1)
  Nz = size(w,2)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
  dXdxI = Metric.dXdxI
  DX = Fe.DX
  DZT = Fe.DZT

  tempZ = Cache.Column
  @views DtempZ = Cache.Block[:,:,1]
  @views tempX = Cache.Block[:,:,2]
  temp = tempX
  @views DtempX = Cache.Block[:,:,3]
  @views FluxZ = Cache.Block[:,1,4]
  OPz = OrdPolyZ + 1
  for ix = 1 : Nx 
    for iz = 1 : Nz 
      @views @. tempZ[:,:,iz] = dXdxI[ix,iz,:,:,2,2] * uLGL[ix,iz,:,:] - dXdxI[ix,iz,:,:,2,1] * w[ix,iz,:,:]
    end  
    for iz = 1 : Nz 
      # Momentum  
      @views @. tempX = (dXdxI[ix,iz,:,:,1,2] * uLGL[ix,iz,:,:] - dXdxI[ix,iz,:,:,1,1] * w[ix,iz,:,:])  
      if iz > 1
        @views @. FluxZ = 0.5 * (tempZ[:,1,iz] - tempZ[:,OPz,iz-1])
        @views @. FuLGL[ix,iz,:,1] += FluxZ * RhoLGL[ix,iz,:,1] * w[ix,iz,:,1]
        @views @. Fw[ix,iz,:,1] -= FluxZ * RhoLGL[ix,iz,:,1] * uLGL[ix,iz,:,1]
      end
      if iz < Nz
        @views @. FluxZ = 0.5 *(tempZ[:,1,iz+1] - tempZ[:,OPz,iz])
        @views @. FuLGL[ix,iz,:,OPz] += FluxZ * RhoLGL[ix,iz,:,OPz] * w[ix,iz,:,OPz]
        @views @. Fw[ix,iz,:,OPz] -= FluxZ * RhoLGL[ix,iz,:,OPz] * uLGL[ix,iz,:,OPz]
      end    
      mul!(DtempX,DX,tempX)
      @views mul!(DtempZ,tempZ[:,:,iz],DZT)
      @views @. temp = RhoLGL[ix,iz,:,:] * (DtempX + DtempZ)
      @views @. FuLGL[ix,iz,:,:] -= w[ix,iz,:,:] * temp
      @views @. Fw[ix,iz,:,:] += uLGL[ix,iz,:,:] * temp
    end 
  end 
end

