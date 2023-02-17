function DSS!(c,Rho,JC)
  Nz=size(c,2)
  Nx=size(c,1)
  OPx=size(c,3)
  OPz=size(c,4)
  cI=zeros(OPz)
  if OPx > 2
    @views @. c[:,:,2:OPx-1,:] /= (JC[:,:,2:OPx-1,:] * Rho[:,:,2:OPx-1,:]) 
  end
  for ix = 2 : Nx
    for iz = 1 : Nz  
      @views @. cI = (c[ix-1,iz,OPx,:] + c[ix,iz,1,:] ) / 
        (JC[ix-1,iz,OPx,:] * Rho[ix-1,iz,OPx,:] + JC[ix,iz,1,:] * Rho[ix,iz,1,:])
      @views @. c[ix-1,iz,OPx,:] = cI
      @views @. c[ix,iz,1,:] = cI
    end  
  end
  for iz = 1 : Nz  
    @views @. cI = (c[1,iz,1,:] + c[Nx,iz,OPx,:] ) / 
      (JC[1,iz,1,:] * Rho[1,iz,1,:] + JC[Nx,iz,OPx,:] * Rho[Nx,iz,OPx,:])
    @views @. c[1,iz,1,:] = cI
    @views @. c[Nx,iz,OPx,:] = cI
  end  
end

function DSS!(c,JC)
  Nz=size(c,2)
  Nx=size(c,1)
  OPx=size(c,3)
  OPz=size(c,4)
  cI=zeros(OPz)
  if OPx > 2
    @views @. c[:,:,2:OPx-1,:] /= JC[:,:,2:OPx-1,:] 
  end
  for ix = 2 : Nx
    for iz = 1 : Nz  
      @views @. cI = (c[ix-1,iz,OPx,:] + c[ix,iz,1,:] ) / 
        (JC[ix-1,iz,OPx,:] + JC[ix,iz,1,:])
      @views @. c[ix-1,iz,OPx,:] = cI
      @views @. c[ix,iz,1,:] = cI
    end  
  end
  for iz = 1 : Nz  
    @views @. cI = (c[1,iz,1,:] + c[Nx,iz,OPx,:] ) / 
      (JC[1,iz,1,:] + JC[Nx,iz,OPx,:])
    @views @. c[Nx,iz,OPx,:] = cI
    @views @. c[1,iz,1,:] = cI
  end  
end

function Average!(c)
  Nz=size(c,2)
  Nx=size(c,1)
  OPx=size(c,3)
  OPz=size(c,4)
  cIL=zeros(OPz)
  cIR=zeros(OPz)
  cI=zeros(OPz)
  for i=2:Nx
    for j = 1 : Nz  
      @views cIL = c[i-1,j,OPx,:]  
      @views cIR = c[i,j,1,:]  
      @views @. cI = 0.5 * (cIL + cIR) 
      ccI = cI  
      @views @. c[i-1,j,OPx,:] = ccI
      @views @. c[i,j,1,:] = ccI
    end  
  end
  for j = 1 : Nz  
    @views cIL = c[Nx,j,OPx,:]  
    @views cIR = c[1,j,1,:]  
    @views @. cI = 0.5 * (cIL + cIR ) 
    ccI =  cI  
    @views @. c[Nx,j,OPx,:] = ccI
    @views @. c[1,j,1,:] = ccI
  end  
end

