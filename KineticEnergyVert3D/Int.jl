function IntF(c,Rho,J,Fe)

  I = 0.0
  Nx = size(c,1)
  Nz = size(c,2)
  OPx = size(c,3)
  OPz = size(c,4)
  wx = Fe.wX
  wz = Fe.wZ
  W = wx * wz'
  for ix = 1 : Nx
    for iz = 1 : Nz
      @views I = I + sum(W .* c[ix,iz,:,:] .* Rho[ix,iz,:,:] .* J[ix,iz,:,:])
    end
  end
  return I
end

function IntF(c,J,Fe)

  I = 0.0
  Nx = size(c,1)
  Nz = size(c,2)
  OPx = size(c,3)
  OPz = size(c,4)
  wx = Fe.wX
  wz = Fe.wZ
  W = wx * wz'
  for ix = 1 : Nx
    for iz = 1 : Nz
      @views I = I + sum(W .* c[ix,iz,:,:] .* J[ix,iz,:,:])
    end
  end
  return I
end

function IntF(c,Fe)

  I = 0.0
  Nx = size(c,1)
  Nz = size(c,2)
  OPx = size(c,3)
  OPz = size(c,4)
  wX = Fe.wX
  wZ = Fe.wZ
  W = wX * wZ'
  for ix = 1 : Nx
    for iz = 1 : Nz
      @views I = I + sum(W .* c[ix,iz,:,:])
    end
  end
  return I
end

function IntC(c,Rho,J,Fe)

  Nx = size(c,1)
  Ny = size(c,2)
  Nz = size(c,3)
  OPx = size(c,4)
  OPy = size(c,5)
  OPz = size(c,6)
  wx = Fe.wX
  wy = Fe.wY
  wz = Fe.wZC

  I = 0.0
  for ix = 1 : Nx
    for ix = 1 : Nx
      for iz = 1 : Nz
        for i = 1 : OPx   
          for j = 1 : OPy   
            for k = 1 : OPz   
              I = I + wx[i] * wy[j] * wz[k] * c[ix,iy,iz,i,j,k] * 
                Rho[ix,iy,iz,i,j,k] * J[ix,iy,iz,i,j,k]  
            end
          end
        end
      end
    end
  end
  return I
end

function IntC(c,J,Fe)

  I = 0
  Nx = size(c,1)
  Nz = size(c,2)
  OPx = size(c,3)
  OPz = size(c,4)
  wx = Fe.wX
  wz = Fe.wZC
  W = wx * wz'
  for ix = 1 : Nx
    for iz = 1 : Nz
      @views I = I + sum(W .* c[ix,iz,:,:] .* J[ix,iz,:,:])
    end
  end
  return I
end

