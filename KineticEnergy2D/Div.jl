function DivRho!(FRhoLGL,uLGL,w,RhoLGL,Fe,Metric)
  Nx = size(w,1)
  Nz = size(w,2)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
  dXdxI = Metric.dXdxI
  VolSurfB = Metric.VolSurfB
  VolSurfT = Metric.VolSurfT

  tempX = zeros(OrdPolyX+1,OrdPolyZ+1)
  DtempX = zeros(OrdPolyX+1,OrdPolyZ+1)
  tempZ = zeros(Nz,OrdPolyX+1,OrdPolyZ+1)
  DtempZ = zeros(OrdPolyX+1,OrdPolyZ+1)
  temp = zeros(OrdPolyX+1,OrdPolyZ+1)
  FluxZ = zeros(OrdPolyX+1)
  OPz = OrdPolyZ + 1
  for ix = 1 : Nx 
    for iz = 1 :Nz  
      @views @. tempZ[iz,:,:] = RhoLGL[ix,iz,:,:] * (dXdxI[ix,iz,:,:,2,1] * uLGL[ix,iz,:,:] +
        dXdxI[ix,iz,:,:,2,2] * w[ix,iz,:,:])
    end    
    for iz = 1 :Nz  
      @views @. tempX = RhoLGL[ix,iz,:,:] .* (dXdxI[ix,iz,:,:,1,1] .* uLGL[ix,iz,:,:] +
        dXdxI[ix,iz,:,:,1,2] .* w[ix,iz,:,:])
      mul!(DtempX,Fe.DX,tempX)
      @views mul!(DtempZ,tempZ[iz,:,:],Fe.DZT)
      @views @. FRhoLGL[ix,iz,:,:] -= (DtempX + DtempZ) 
      # Density
      if iz > 1
        @views @. FluxZ = 0.5 * (tempZ[iz,:,1] -tempZ[iz-1,:,OPz]) 
        @views @. FRhoLGL[ix,iz,:,1] -= FluxZ
      end  
      if iz < Nz
        @views @. FluxZ = 0.5 * (tempZ[iz+1,:,1] -tempZ[iz,:,OPz]) 
        @views @. FRhoLGL[ix,iz,:,OPz] -= FluxZ
      end  
    end 
  end 
end

function DivRhoTheta!(FRhoThetaLGL,uLGL,w,RhoThetaLGL,Fe,Metric)
  Nx = size(w,1)
  Nz = size(w,2)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
  dXdxI = Metric.dXdxI
  VolSurfB = Metric.VolSurfB
  VolSurfT = Metric.VolSurfT

  tempX = zeros(OrdPolyX+1,OrdPolyZ+1)
  DtempX = zeros(OrdPolyX+1,OrdPolyZ+1)
  tempZ = zeros(Nz,OrdPolyX+1,OrdPolyZ+1)
  DtempZ = zeros(OrdPolyX+1,OrdPolyZ+1)
  temp = zeros(OrdPolyX+1,OrdPolyZ+1)
  FluxZ = zeros(OrdPolyX+1)
  OPz = OrdPolyZ + 1
  for ix = 1 : Nx 
    for iz = 1 :Nz  
      @views @. tempZ[iz,:,:] = RhoThetaLGL[ix,iz,:,:] * (dXdxI[ix,iz,:,:,2,1] * uLGL[ix,iz,:,:] +
        dXdxI[ix,iz,:,:,2,2] * w[ix,iz,:,:])
    end    
    for iz = 1 :Nz  
      @views @. tempX = RhoThetaLGL[ix,iz,:,:] .* (dXdxI[ix,iz,:,:,1,1] .* uLGL[ix,iz,:,:] +
        dXdxI[ix,iz,:,:,1,2] .* w[ix,iz,:,:])
      mul!(DtempX,Fe.DX,tempX)
      @views mul!(DtempZ,tempZ[iz,:,:],Fe.DZT)
      @views @. FRhoThetaLGL[ix,iz,:,:] -= (DtempX + DtempZ) 
      # Density
      if iz > 1
        @views @. FluxZ = 0.5 * (tempZ[iz,:,1] -tempZ[iz-1,:,OPz]) 
        @views @. FRhoThetaLGL[ix,iz,:,1] -= FluxZ
      end  
      if iz < Nz
        @views @. FluxZ = 0.5 * (tempZ[iz+1,:,1] -tempZ[iz,:,OPz]) 
        @views @. FRhoThetaLGL[ix,iz,:,OPz] -= FluxZ
      end  
    end 
  end 
end

