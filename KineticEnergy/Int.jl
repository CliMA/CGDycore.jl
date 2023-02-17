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
      @views I = I + sum(W .* c[ix,iz,:,:] .* Rho[ix,iz,:,:] .* J[ix,iz,:,:])
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

