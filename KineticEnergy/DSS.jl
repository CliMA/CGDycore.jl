function DSS!(c,J,I12,I21)
  Nz=size(c,2)
  Nx=size(c,1)
  OPx=size(c,3)
  OPz=size(c,4)
  cIL=zeros(OPz+1)
  cIR=zeros(OPz+1)
  cI=zeros(OPz+1)
  for i=2:Nx
    for j = 1 : Nz  
      @views cIL = I12 * c[i-1,j,OPx,:]  
      @views cIR =I12 * c[i,j,1,:]  
      @views @. cI = (cIL + cIR ) / (J[i-1,j,OPx,:] + J[i,j,1,:])
      ccI = I21 * cI  
      @views @. c[i-1,j,OPx,:] = ccI
      @views @. c[i,j,1,:] = ccI
    end  
  end
  for j = 1 : Nz  
    @views cIL = I12 * c[Nx,j,OPx,:]  
    @views cIR =I12 * c[1,j,1,:]  
    @views @. cI = (cIL  + cIR) / (J[Nx,j,OPx,:] + J[1,j,1,:])
    ccI = I21 * cI  
    @views @. c[Nx,j,OPx,:] = ccI
    @views @. c[1,j,1,:] = ccI
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

