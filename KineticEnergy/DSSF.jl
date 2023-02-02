function DSSF!(c,J)
  Nz=size(c,2)
  Nx=size(c,1)
  OPx=size(c,3)
  OPz=size(c,4)
  if OPz > 2
    cL = zeros(OPz-2)
    cR = zeros(OPz-2)
    cA = zeros(OPz-2)
    for j=2:Nx
      for i=1:Nz
        @views @. cL = c[j-1,i,OPx,2:OPz-1]  
        @views @. cR = c[j,i,1,2:OPz-1]  
        @views @. cA = (cL + cR) / (J[j-1,i,OPx,2:OPz-1] + J[j,i,1,2:OPz-1])
        @views @. c[j-1,i,OPx,2:OPz-1] = cA
        @views @. c[j,i,1,2:OPz-1] = cA
      end
    end
    for i=1:Nz
      @views @. cL = c[Nx,i,OPx,2:OPz-1]  
      @views @. cR = c[1,i,1,2:OPz-1]  
      @views @. cA = (cL + cR) / (J[Nx,i,OPx,2:OPz-1] + J[1,i,1,2:OPz-1])
      @views @. c[Nx,i,OPx,2:OPz-1] = cA
      @views @. c[1,i,1,2:OPz-1] = cA
    end
  end
  if OPx > 2
    cL = zeros(OPx-2)
    cR = zeros(OPx-2)
    cA = zeros(OPx-2)
    for j=1:Nx
      for i=2:Nz
        @views @. cL = c[j,i-1,2:OPx-1,OPz]  
        @views @. cR = c[j,i,2:OPx-1,1]  
        @views @. cA = (cL + cR) / (J[j,i-1,2:OPx-1,OPz] + J[j,i,2:OPx-1,1])
        @views @. c[j,i-1,2:OPx-1,OPz] = cA
        @views @. c[j,i,2:OPx-1,1] = cA
      end
    end
  end
  for j=2:Nx
    for i=2:Nz
      cv = (c[j-1,i-1,OPx,OPz]  + c[j-1,i,OPx,1] +  
          + c[j,i-1,1,OPz]  + c[j,i,1,1]) /   
        (J[j-1,i-1,OPx,OPz] + J[j-1,i,OPx,1] +
         J[j,i-1,1,OPz] + J[j,i,1,1])  
      c[j-1,i-1,OPx,OPz] = cv
      c[j-1,i,OPx,1] = cv
      c[j,i-1,1,OPz] = cv
      c[j,i,1,1] = cv
    end
    cv = (c[j-1,1,OPx,1] + c[j,1,1,1]) /   
         (J[j-1,1,OPx,1] + J[j,1,1,1])   
    c[j-1,1,OPx,1] = cv
    c[j,1,1,1] = cv
    cv = (c[j-1,Nz,OPx,OPz] + c[j,Nz,1,OPz]) /   
         (J[j-1,Nz,OPx,OPz] + J[j,Nz,1,OPz])   
    c[j-1,Nz,OPx,OPz] = cv
    c[j,Nz,1,OPz] = cv
  end
  for i=2:Nz
    cv = (c[Nx,i-1,OPx,OPz]  + c[Nx,i,OPx,1] +
        + c[1,i-1,1,OPz]  + c[1,i,1,1]) /   
      (J[Nx,i-1,OPx,OPz] + J[Nx,i,OPx,1] +
       J[1,i-1,1,OPz] + J[1,i,1,1])  
    c[Nx,i-1,OPx,OPz] = cv
    c[Nx,i,OPx,1] = cv
    c[1,i-1,1,OPz] = cv
    c[1,i,1,1] = cv
  end
  cv = (c[Nx,1,OPx,1] + c[1,1,1,1]) /   
       (J[Nx,1,OPx,1] + J[1,1,1,1])   
  c[Nx,1,OPx,1] = cv
  c[1,1,1,1] = cv
  cv = (c[Nx,Nz,OPx,OPz] + c[1,Nz,1,OPz]) /   
       (J[Nx,Nz,OPx,OPz] + J[1,Nz,1,OPz])   
  c[Nx,Nz,OPx,OPz] = cv
  c[1,Nz,1,OPz] = cv
end  

function AverageF!(c,J)
  Nz=size(c,2)
  Nx=size(c,1)
  OPx=size(c,3)
  OPz=size(c,4)
  if OPz > 2
    cL = zeros(OPz-2)
    cR = zeros(OPz-2)
    cA = zeros(OPz-2)
    for j=2:Nx
      for i=1:Nz
        cL = c[j-1,i,OPx,2:OPz-1]  
        cR = c[j,i,1,2:OPz-1]  
        cA = 0.5 * (cL + cR)
        c[j-1,i,OPx,2:OPz-1] = cA
        c[j,i,1,2:OPz-1] = cA
      end
    end
    for i=1:Nz
      cL = c[Nx,i,OPx,2:OPz-1]  
      cR = c[1,i,1,2:OPz-1]  
      cA = 0.5 * (cL + cR)
      c[Nx,i,OPx,2:OPz-1] = cA
      c[1,i,1,2:OPz-1] = cA
    end
  end
  if OPx > 2
    cL = zeros(OPx-2)
    cR = zeros(OPx-2)
    cA = zeros(OPx-2)
    for j=1:Nx
      for i=2:Nz
        cL = c[j,i,2:OPx-1,OPz]  
        cR = c[j,i+1,2:OPx-1,1]  
        cA = 0.5 * (cL + cR) 
        c[j,i-1,2:OPx-1,OPz] = cA
        c[j,i,2:OPx-1,1] = cA
      end
    end
  end
  for j=2:Nx
    for i=2:Nz
      cv = 0.25 * (c[j-1,i-1,OPx,OPz]  + c[j-1,i,OPx,1]  
          + c[j,i-1,1,OPz]  + c[j,i,1,1])    
      c[j-1,i-1,OPx,OPz] = cv
      c[j-1,i,OPx,1] = cv
      c[j,i-1,1,OPz] = cv
      c[j,i,1,1] = cv
    end
    cv = 0.5 * (c[j-1,1,OPx,1] + c[j,1,1,1])
    c[j-1,1,OPx,1] = cv
    c[j,1,1,1] = cv
    cv = 0.5 * (c[j-1,Nz,OPx,OPz] + c[j,Nz,1,OPz])
    c[j-1,Nz,OPx,OPz] = cv
    c[j,Nz,1,OPz] = cv

  end
  for i=2:Nz
    cv = 0.25 * (c[Nx,i-1,OPx,OPz]  + c[Nx,i,OPx,1]  
        + c[1,i-1,1,OPz]  + c[1,i,1,1])    
    c[Nx,i-1,OPx,OPz] = cv
    c[Nx,i,OPx,1] = cv
    c[1,i-1,1,OPz] = cv
    c[1,i,1,1] = cv
  end
  cv = 0.5 * (c[Nx,1,OPx,1] + c[1,1,1,1])    
  c[Nx,1,OPx,1] = cv
  c[1,1,1,1] = cv
  cv = 0.5 * (c[Nx,Nz,OPx,OPz] + c[1,Nz,1,OPz]) 
  c[Nx,Nz,OPx,OPz] = cv
  c[1,Nz,1,OPz] = cv
end  


