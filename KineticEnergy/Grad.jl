function RhoGrad!(FuLGL,Fw,pLGL,RhoLGL,Fe,Metric)
  Nx = size(pLGL,1)
  Nz = size(pLGL,2)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
  dXdxI = Metric.dXdxI

  DtempX = zeros(OrdPolyX+1,OrdPolyZ+1)
  DtempZ = zeros(OrdPolyX+1,OrdPolyZ+1)
  FluxZ = zeros(OrdPolyX+1)
  GradZ = zeros(OrdPolyX+1)
  OPz = OrdPolyZ + 1

  for ix = 1 : Nx 
    for iz = 1 : Nz 
      @views mul!(DtempX,Fe.DX,pLGL[ix,iz,:,:])
      @views mul!(DtempZ,pLGL[ix,iz,:,:],Fe.DZT)
      @views @. FuLGL[ix,iz,:,:] -= RhoLGL[ix,iz,:,:] * 
        (dXdxI[ix,iz,:,:,1,1] * DtempX + dXdxI[ix,iz,:,:,2,1] * DtempZ)
      @views @. Fw[ix,iz,:,:] -= RhoLGL[ix,iz,:,:] * 
        (dXdxI[ix,iz,:,:,1,2] * DtempX + dXdxI[ix,iz,:,:,2,2] * DtempZ)
      if iz > 1
        @views @. GradZ = 0.5 * (pLGL[ix,iz,:,1] - pLGL[ix,iz-1,:,OPz])  
        @views @. FluxZ = GradZ * dXdxI[ix,iz,:,1,2,1] * RhoLGL[ix,iz,:,1] 
        @views @. FuLGL[ix,iz,:,1] -= FluxZ 
        @views @. FluxZ = GradZ * dXdxI[ix,iz,:,1,2,2] * RhoLGL[ix,iz,:,1] 
        @views @. Fw[ix,iz,:,1] -= FluxZ 
      end    
      if iz < Nz 
        @views @. GradZ = 0.5 * (pLGL[ix,iz+1,:,1] - pLGL[ix,iz,:,OPz])  
        @views @. FluxZ = GradZ * dXdxI[ix,iz,:,OPz,2,1] * RhoLGL[ix,iz,:,OPz]
        @views @. FuLGL[ix,iz,:,OPz] -= FluxZ 
        @views @. FluxZ =  GradZ * dXdxI[ix,iz,:,OPz,2,2] * RhoLGL[ix,iz,:,OPz]
        @views @. Fw[ix,iz,:,OPz] -= FluxZ 
      end    
    end 
  end 
end

function Grad!(FuLGL,Fw,pLGL,Fe,Metric)
  Nx = size(pLGL,1)
  Nz = size(pLGL,2)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
  dXdxI = Metric.dXdxI

  DtempX = zeros(OrdPolyX+1,OrdPolyZ+1)
  DtempZ = zeros(OrdPolyX+1,OrdPolyZ+1)
  FluxZ = zeros(OrdPolyX+1)
  GradZ = zeros(OrdPolyX+1)
  OPz = OrdPolyZ + 1

  for ix = 1 : Nx 
    for iz = 1 : Nz 
      @views mul!(DtempX,Fe.DX,pLGL[ix,iz,:,:])
      @views mul!(DtempZ,pLGL[ix,iz,:,:],Fe.DZT)
      @views @. FuLGL[ix,iz,:,:] -= (dXdxI[ix,iz,:,:,1,1] * DtempX + dXdxI[ix,iz,:,:,2,1] * DtempZ)
      @views @. Fw[ix,iz,:,:] -= (dXdxI[ix,iz,:,:,1,2] * DtempX + dXdxI[ix,iz,:,:,2,2] * DtempZ)
      if iz > 1
        @views @. GradZ = 0.5 * (pLGL[ix,iz,:,1] - pLGL[ix,iz-1,:,OPz])  
        @views @. FluxZ = GradZ * dXdxI[ix,iz,:,1,2,1]
        @views @. FuLGL[ix,iz,:,1] -= FluxZ 
        @views @. FluxZ = GradZ * dXdxI[ix,iz,:,1,2,2]
        @views @. Fw[ix,iz,:,1] -= FluxZ 
      end    
      if iz < Nz 
        @views @. GradZ = 0.5 * (pLGL[ix,iz+1,:,1] - pLGL[ix,iz,:,OPz])  
        @views @. FluxZ = GradZ * dXdxI[ix,iz,:,OPz,2,1]
        @views @. FuLGL[ix,iz,:,OPz] -= FluxZ 
        @views @. FluxZ =  GradZ * dXdxI[ix,iz,:,OPz,2,2]
        @views @. Fw[ix,iz,:,OPz] -= FluxZ 
      end    
    end 
  end 
end

