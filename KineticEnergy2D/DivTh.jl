function Div!(FThLGL,ThLGL,uLGL,w,RhoLGL,Fe,Metric)
  Nx = size(w,1)
  Nz = size(w,2)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
  dXdxI = Metric.dXdxI
  VolSurfB = Metric.VolSurfB
  VolSurfT = Metric.VolSurfT

  FRhoLGL .= 0.0

  tempX = zeros(OrdPolyX+1,OrdPolyZ+1)
  DtempX = zeros(OrdPolyX+1,OrdPolyZ+1)
  tempZ = zeros(Nz,OrdPolyX+1,OrdPolyZ+1)
  DtempZ = zeros(OrdPolyX+1,OrdPolyZ+1)
  temp = zeros(OrdPolyX+1,OrdPolyZ+1)
  FluxZ = zeros(OrdPolyX+1)
  OPz = OrdPolyZ + 1
  for ix = 1 : Nx 
    for iz = 1 :Nz  
      @views @. tempZ[iz,:,:] = ThLGL[ix,iz,:,:] * (dXdxI[ix,iz,:,:,2,1] * uLGL[ix,iz,:,:] +
        dXdxI[ix,iz,:,:,2,2] * w[ix,iz,:,:])
    end    
    for iz = 1 :Nz  
      # Density
      @views @. tempX = ThLGL[ix,iz,:,:] .* (dXdxI[ix,iz,:,:,1,1] .* uLGL[ix,iz,:,:] +
        dXdxI[ix,iz,:,:,1,2] .* w[ix,iz,:,:])
      if iz > 1
        @views @. FluxZ = (temp[iz,:,1] -temp[iz-1,:,OPz]) * VolSurfB[ix,iz,:]
        @views @. FThLGL[ix,iz-1,:,OPz] -= FluxZ
        @views @. FThLGL[ix,iz,:,1] += FluxZ
      end  
      mul!(DtempX,Fe.DX,tempX)
      @views mul!(DtempZ,tempZ[iz,:,:],Fe.DZ')
      @views @. FThLGL[ix,iz,:,:] -= (tempX + tempZ) 
    end 
  end 
end

