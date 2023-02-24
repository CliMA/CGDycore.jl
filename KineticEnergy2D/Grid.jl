mutable struct Grid
  Nx::Int
  Nz::Int
  Lx::Float64
  Lz::Float64
  x0::Float64
  xP::Array{Float64, 4}
  zP::Array{Float64, 4}
end

function Grid(Nx,Nz,x0,Lx,Lz,Oro,Fe,Param)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
  #Grid 
  xP=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  @views @. xP[1,:,1,:] = x0
  dx=Lx/Nx
  for ix = 1 : Nx
    for iz = 1 : Nz  
      for i = 2 : OrdPolyX + 1
        for j = 1 : OrdPolyZ + 1
          xP[ix,iz,i,j] = xP[ix,iz,1,j] + (1 + Fe.xe[i]) / 2 * dx
        end
      end
    end
    if ix < Nx
      @views @. xP[ix+1,:,1,:] = xP[ix,:,OrdPolyX+1,:]
    end
  end  

  @views @. xP[Nx,:,OrdPolyX+1,:] = Lx + x0

  zP=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  for ix = 1 : Nx
    for i = 1 : OrdPolyX + 1  
      zP[ix,1,i,1]=Oro(xP[ix,1,i,1],Param)
      dzLoc=(Lz - zP[ix,1,i,1]) / Nz
      for iz = 1 : Nz
        for j = 1 : OrdPolyZ + 1  
          zP[ix,iz,i,j]=zP[ix,iz,i,1]+(1+Fe.ze[j])/2*dzLoc   
        end  
        if iz < Nz
          zP[ix,iz+1,i,1] = zP[ix,iz,i,OrdPolyZ+1]  
        end  
      end  
      zP[ix,Nz,i,OrdPolyZ+1] = Lz
    end
  end
  return Grid(
    Nx,
    Nz,
    Lx,
    Lz,
    x0,
    xP,
    zP,
  )
end  
