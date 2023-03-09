mutable struct Grid
  Nx::Int
  Ny::Int
  Nz::Int
  Lx::Float64
  Ly::Float64
  Lz::Float64
  x0::Float64
  y0::Float64
  xP::Array{Float64, 6}
  yP::Array{Float64, 6}
  zP::Array{Float64, 6}
end

function Grid(Nx,Ny,Nz,x0,y0,Lx,Ly,Lz,Oro,Fe,Param)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  #Grid 
  xP = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1)
  @views @. xP[1,:,:,1,:,:] = x0
  dx = Lx / Nx
  for ix = 1 : Nx
    for iy = 1 : Ny
      for iz = 1 : Nz  
        for i = 2 : OrdPolyX + 1
          for j = 1 : OrdPolyY + 1
            for k = 1 : OrdPolyZ + 1
              xP[ix,iy,iz,i,j,k] = xP[ix,iy,iz,1,j,k] + (1 + Fe.xe[i]) / 2 * dx
            end
          end
        end
      end
    end
    if ix < Nx
      @views @. xP[ix+1,:,:,1,:,:] = xP[ix,:,:,OrdPolyX+1,:,:]
    end
  end  
  @views @. xP[Nx,:,:,OrdPolyX+1,:,:] = Lx + x0

  yP = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1)
  @views @. yP[:,1,:,:,1,:] = y0
  dy = Ly / Ny
  for ix = 1 : Nx
    for iy = 1 : Ny
      for iz = 1 : Nz  
        for i = 1 : OrdPolyX + 1
          for j = 2 : OrdPolyY + 1
            for k = 1 : OrdPolyZ + 1
              yP[ix,iy,iz,i,j,k] = yP[ix,iy,iz,i,1,k] + (1 + Fe.ye[j]) / 2 * dy
            end
          end
        end
      end
      if iy < Ny
        @views @. yP[:,iy+1,:,:,1,:] = yP[:,iy,:,:,OrdPolyY+1,:]
      end
    end
  end  
  @views @. yP[:,Ny,:,:,OrdPolyY+1,:] = Ly + y0

  zP=zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1)
  for ix = 1 : Nx
    for iy = 1 : Ny
      for i = 1 : OrdPolyX + 1  
        for j = 1 : OrdPolyY + 1  
          zP[ix,iy,1,i,j,1]=Oro(xP[ix,iy,1,i,j,1],yP[ix,iy,1,i,j,1],Param)
          dzLoc=(Lz - zP[ix,iy,1,i,j,1]) / Nz
          for iz = 1 : Nz
            for k = 1 : OrdPolyZ + 1  
              zP[ix,iy,iz,i,j,k]=zP[ix,iy,iz,i,j,1]+(1+Fe.ze[k])/2*dzLoc 
            end  
            if iz < Nz
              zP[ix,iy,iz+1,i,j,1] = zP[ix,iy,iz,i,j,OrdPolyZ+1]  
            end  
          end  
          zP[ix,iy,Nz,i,j,OrdPolyZ+1] = Lz
        end  
      end  
    end
  end
  return Grid(
    Nx,
    Ny,
    Nz,
    Lx,
    Ly,
    Lz,
    x0,
    y0,
    xP,
    yP,
    zP,
  )
end  
